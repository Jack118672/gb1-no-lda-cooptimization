from __future__ import annotations

import json
import os
import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd
import torch


ROOT = Path.cwd()
FITNESS_PATH = ROOT / "gb1_fitness_singles.csv"
OUT_PATH = ROOT / "gb1_lm_scores.csv"
CACHE_HUB = ROOT / "esm_cache" / "hub"

# Olson/MAVE GB1 table uses positions 2-56 and omits the N-terminal M.
WT_SEQUENCE_BINDING_TABLE = "QYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
WT_POSITION_START = 2
MODEL_NAME = "esm2_t6_8M_UR50D"


def mutation_label(pos: int, wt: str, mut: str) -> str:
    return f"{wt}{int(pos)}{mut}"


def main() -> None:
    if not FITNESS_PATH.exists():
        raise FileNotFoundError(FITNESS_PATH)

    # ESM uses torch.hub's cache layout. The model is downloaded/cached here.
    if CACHE_HUB.exists():
        torch.hub.set_dir(str(CACHE_HUB))
        os.environ["TORCH_HOME"] = str(ROOT / "esm_cache")

    # Recent PyTorch builds default to a stricter checkpoint loader. The cached
    # ESM1b checkpoint stores argparse.Namespace metadata from the trusted Meta
    # ESM release, so we allow-list that class before esm.pretrained loads it.
    if hasattr(torch.serialization, "add_safe_globals"):
        torch.serialization.add_safe_globals([argparse.Namespace])

    import esm

    singles = pd.read_csv(FITNESS_PATH)
    singles["mutation"] = [
        mutation_label(pos, wt, mut)
        for pos, wt, mut in singles[["pos", "wt", "mut"]].itertuples(index=False, name=None)
    ]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, alphabet = esm.pretrained.load_model_and_alphabet(MODEL_NAME)
    model.eval()
    model = model.to(device)

    batch_converter = alphabet.get_batch_converter()
    _, _, batch_tokens = batch_converter([("GB1_binding_table_WT", WT_SEQUENCE_BINDING_TABLE)])
    batch_tokens = batch_tokens.to(device)

    with torch.no_grad():
        logits = model(batch_tokens)["logits"]
        token_log_probs = torch.log_softmax(logits, dim=-1).cpu()[0]

    rows = []
    for row in singles.itertuples(index=False):
        pos = int(row.pos)
        wt = row.wt
        mut = row.mut
        seq_index = pos - WT_POSITION_START
        if seq_index < 0 or seq_index >= len(WT_SEQUENCE_BINDING_TABLE):
            continue
        sequence_wt = WT_SEQUENCE_BINDING_TABLE[seq_index]
        if sequence_wt != wt:
            raise ValueError(
                f"WT mismatch for {mutation_label(pos, wt, mut)}: "
                f"sequence has {sequence_wt} at dataset position {pos}"
            )

        # Add 1 because ESM token position 0 is BOS.
        token_position = seq_index + 1
        wt_token = alphabet.get_idx(wt)
        mut_token = alphabet.get_idx(mut)
        raw_llr = (
            token_log_probs[token_position, mut_token]
            - token_log_probs[token_position, wt_token]
        ).item()
        rows.append({
            "mutation": mutation_label(pos, wt, mut),
            "pos": pos,
            "wt": wt,
            "mut": mut,
            "esm_model": MODEL_NAME,
            "esm_wt_marginal_llr": raw_llr,
        })

    scores = pd.DataFrame(rows)
    mean = scores["esm_wt_marginal_llr"].mean()
    std = scores["esm_wt_marginal_llr"].std(ddof=0)
    scores["lm_score"] = (scores["esm_wt_marginal_llr"] - mean) / std
    scores = scores.sort_values("lm_score", ascending=False)
    scores.to_csv(OUT_PATH, index=False)

    metadata = {
        "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "output": str(OUT_PATH),
        "model": MODEL_NAME,
        "device": str(device),
        "wt_sequence_binding_table": WT_SEQUENCE_BINDING_TABLE,
        "wt_position_start": WT_POSITION_START,
        "scoring": "ESM WT-marginal log-likelihood ratio: log P(mut | WT context) - log P(wt | WT context)",
        "lm_score": "z-scored esm_wt_marginal_llr; higher means more ESM-plausible",
        "rows": int(len(scores)),
    }
    (ROOT / "gb1_lm_scores_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
