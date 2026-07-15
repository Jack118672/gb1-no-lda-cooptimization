from __future__ import annotations

import itertools
import json
import re
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path.cwd()
OUT = ROOT / "GB1_NO_LDA_RESULTS"
FIG_DIR = OUT / "figures"
INPUT_DIR = OUT / "inputs"
TABLE_DIR = OUT / "tables"

SINGLE_FITNESS_PATH = ROOT / "gb1_fitness_singles.csv"
SINGLE_STABILITY_PATH = ROOT / "gb1_stability_deltaG.csv"
DOUBLE_FITNESS_PATH = ROOT / "gb1_fitness_doubles.csv"
STABILITY_SOURCE_PATH = ROOT / "Data_tables_for_figs" / "dG_extdG_data_Fig1.csv"
OPTIONAL_REAL_LM_SCORES = ROOT / "gb1_lm_scores.csv"

WT_SEQUENCE_BINDING_TABLE = "QYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
WT_BINDING_FITNESS = 1.0
TARGET_ROUND1_CANDIDATES = 20
MAX_ROUND2_DOUBLES_TO_HIGHLIGHT = 30

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

HYDROPATHY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

VOLUME = {
    "A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9,
    "G": 60.1, "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7,
    "M": 162.9, "N": 114.1, "P": 112.7, "Q": 143.8, "R": 173.4,
    "S": 89.0, "T": 116.1, "V": 140.0, "W": 227.8, "Y": 193.6,
}

CHARGE = {
    "D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1,
    **{aa: 0.0 for aa in AA_ORDER if aa not in "DEKRH"},
}

BLOSUM62_TEXT = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""


def parse_blosum62(text: str) -> dict[tuple[str, str], float]:
    lines = [line.strip() for line in text.strip().splitlines()]
    headers = lines[0].split()
    scores = {}
    for line in lines[1:]:
        parts = line.split()
        row = parts[0]
        for col, value in zip(headers, parts[1:]):
            scores[(row, col)] = float(value)
    return scores


BLOSUM62 = parse_blosum62(BLOSUM62_TEXT)


def mutation_label(pos: int, wt: str, mut: str) -> str:
    return f"{wt}{int(pos)}{mut}"


def lm_plausibility_proxy(row: pd.Series) -> float:
    wt = row["wt"]
    mut = row["mut"]
    substitution_similarity = BLOSUM62[(wt, mut)] / 11.0
    hydropathy_similarity = 1.0 - abs(HYDROPATHY[mut] - HYDROPATHY[wt]) / 9.0
    volume_similarity = 1.0 - abs(VOLUME[mut] - VOLUME[wt]) / 200.0
    charge_similarity = 1.0 - abs(CHARGE[mut] - CHARGE[wt]) / 2.0
    return (
        0.50 * substitution_similarity
        + 0.25 * hydropathy_similarity
        + 0.15 * volume_similarity
        + 0.10 * charge_similarity
    )


def pareto_mask_2d(df: pd.DataFrame, x_col: str, y_col: str) -> np.ndarray:
    """Efficient Pareto front for two maximize objectives."""
    work = df[[x_col, y_col]].copy()
    work["_original_index"] = np.arange(len(work))
    work = work.sort_values([x_col, y_col], ascending=[False, False])

    keep = np.zeros(len(df), dtype=bool)
    best_y_from_higher_x = -np.inf
    for _, group in work.groupby(x_col, sort=False):
        group_max_y = group[y_col].max()
        if group_max_y > best_y_from_higher_x:
            group_keep = group[group[y_col] == group_max_y]
            keep[group_keep["_original_index"].to_numpy()] = True
        best_y_from_higher_x = max(best_y_from_higher_x, group_max_y)
    return keep


def ensure_clean_dirs() -> None:
    if OUT.exists():
        for path in sorted(OUT.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir() and path not in {FIG_DIR, INPUT_DIR, TABLE_DIR}:
                try:
                    path.rmdir()
                except OSError:
                    pass
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)


def crosscheck_stability_metric(stability: pd.DataFrame) -> pd.DataFrame:
    if not STABILITY_SOURCE_PATH.exists():
        return pd.DataFrame([{
            "check": "source file exists",
            "value": "missing",
            "interpretation": "Could not cross-check against dG_extdG_data_Fig1.csv.",
        }])

    source = pd.read_csv(STABILITY_SOURCE_PATH)
    source = source[source["PDB_ID"].eq("1PGA")].copy()
    parsed_rows = []
    for _, row in source.iterrows():
        match = re.search(r"_([A-Z])(\d+)([A-Z])$", str(row["name"]))
        if not match:
            continue
        wt, pos, mut = match.groups()
        parsed_rows.append({
            "pos": int(pos),
            "wt": wt,
            "mut": mut,
            "source_deltaG": row["deltaG"],
            "source_prev_dG": row["prev_dG"],
        })
    parsed = pd.DataFrame(parsed_rows)
    compare = stability.merge(parsed, on=["pos", "wt", "mut"], how="inner")
    compare["abs_diff"] = (compare["stability_deltaG"] - compare["source_deltaG"]).abs()

    return pd.DataFrame([
        {
            "check": "local stability column source",
            "value": "gb1_stability_deltaG.csv column ddg equals source deltaG",
            "interpretation": "Use as the source-reported stability scale; larger values are treated as more stable in this dataset.",
        },
        {
            "check": "matched source rows",
            "value": int(len(compare)),
            "interpretation": "Rows matched by pos/wt/mut against source figure table.",
        },
        {
            "check": "max absolute difference vs source deltaG",
            "value": float(compare["abs_diff"].max()),
            "interpretation": "0 means the extracted values are copied directly from source deltaG.",
        },
    ])


def load_and_score_singles() -> tuple[pd.DataFrame, pd.DataFrame, str]:
    fitness = pd.read_csv(SINGLE_FITNESS_PATH)
    stability = pd.read_csv(SINGLE_STABILITY_PATH).rename(columns={"ddg": "stability_deltaG"})

    fitness["mutation"] = [
        mutation_label(pos, wt, mut)
        for pos, wt, mut in fitness[["pos", "wt", "mut"]].itertuples(index=False, name=None)
    ]
    stability["mutation"] = [
        mutation_label(pos, wt, mut)
        for pos, wt, mut in stability[["pos", "wt", "mut"]].itertuples(index=False, name=None)
    ]
    checks = crosscheck_stability_metric(stability)

    singles = fitness.merge(
        stability[["pos", "wt", "mut", "mutation", "stability_deltaG"]],
        on=["pos", "wt", "mut", "mutation"],
        how="inner",
        validate="one_to_one",
    )

    if OPTIONAL_REAL_LM_SCORES.exists():
        lm = pd.read_csv(OPTIONAL_REAL_LM_SCORES)
        singles = singles.merge(lm[["mutation", "lm_score"]], on="mutation", how="left")
        if singles["lm_score"].isna().any():
            raise ValueError("gb1_lm_scores.csv did not cover every labeled single mutant.")
        lm_score_source = "gb1_lm_scores.csv"
    else:
        raw = singles.apply(lm_plausibility_proxy, axis=1)
        singles["lm_score"] = (raw - raw.mean()) / raw.std(ddof=0)
        lm_score_source = "BLOSUM62 + hydropathy/volume/charge similarity proxy"

    singles["measured_pareto_front"] = pareto_mask_2d(singles, "fitness", "stability_deltaG")
    return singles, checks, lm_score_source


def select_round1_candidates(singles: pd.DataFrame) -> pd.DataFrame:
    threshold = singles["lm_score"].nlargest(TARGET_ROUND1_CANDIDATES).min()
    selected = singles[singles["lm_score"] >= threshold].copy()
    selected = selected.sort_values(["lm_score", "fitness"], ascending=[False, False]).reset_index(drop=True)
    selected["binding_improves_vs_wt"] = selected["fitness"] > WT_BINDING_FITNESS
    selected["on_measured_pareto_front"] = selected["measured_pareto_front"]
    return selected


def build_round2_candidates(selected: pd.DataFrame, doubles: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pairs = []
    for first, second in itertools.combinations(selected.to_dict("records"), 2):
        if int(first["pos"]) == int(second["pos"]):
            continue
        a, b = sorted([first, second], key=lambda row: int(row["pos"]))
        pairs.append({
            "double_mutation": f"{a['mutation']}+{b['mutation']}",
            "pos1": int(a["pos"]),
            "mut1": a["mut"],
            "mutation_1": a["mutation"],
            "pos2": int(b["pos"]),
            "mut2": b["mut"],
            "mutation_2": b["mutation"],
            "single1_lm_score": a["lm_score"],
            "single2_lm_score": b["lm_score"],
            "double_lm_score_additive": a["lm_score"] + b["lm_score"],
            "single1_deltaG": a["stability_deltaG"],
            "single2_deltaG": b["stability_deltaG"],
            "double_stability_deltaG_additive": a["stability_deltaG"] + b["stability_deltaG"],
        })

    generated = pd.DataFrame(pairs)
    doubles_clean = doubles.rename(columns={"fitness": "measured_double_fitness"}).copy()
    doubles_clean["pos1"] = doubles_clean["pos1"].astype(int)
    doubles_clean["pos2"] = doubles_clean["pos2"].astype(int)

    generated = generated.merge(
        doubles_clean[["pos1", "mut1", "pos2", "mut2", "measured_double_fitness"]],
        on=["pos1", "mut1", "pos2", "mut2"],
        how="left",
    )
    generated["has_measured_double_fitness"] = generated["measured_double_fitness"].notna()
    generated["binding_improves_vs_wt"] = generated["measured_double_fitness"] > WT_BINDING_FITNESS

    highlighted = (
        generated[generated["has_measured_double_fitness"]]
        .sort_values(["double_lm_score_additive", "measured_double_fitness"], ascending=[False, False])
        .head(MAX_ROUND2_DOUBLES_TO_HIGHLIGHT)
        .copy()
    )
    return generated, highlighted


def build_double_validation_background(singles: pd.DataFrame, doubles: pd.DataFrame) -> pd.DataFrame:
    doubles_clean = doubles.rename(columns={"fitness": "measured_double_fitness"}).copy()
    doubles_clean["pos1"] = doubles_clean["pos1"].astype(int)
    doubles_clean["pos2"] = doubles_clean["pos2"].astype(int)

    single_stability_1 = singles[["pos", "mut", "stability_deltaG"]].rename(
        columns={"pos": "pos1", "mut": "mut1", "stability_deltaG": "single1_deltaG"}
    )
    single_stability_2 = singles[["pos", "mut", "stability_deltaG"]].rename(
        columns={"pos": "pos2", "mut": "mut2", "stability_deltaG": "single2_deltaG"}
    )

    background = doubles_clean.merge(single_stability_1, on=["pos1", "mut1"], how="inner")
    background = background.merge(single_stability_2, on=["pos2", "mut2"], how="inner")
    background["double_stability_deltaG_additive"] = (
        background["single1_deltaG"] + background["single2_deltaG"]
    )
    background["validation_pareto_front"] = pareto_mask_2d(
        background, "measured_double_fitness", "double_stability_deltaG_additive"
    )
    return background


def plot_round1_clean(singles: pd.DataFrame, selected: pd.DataFrame) -> Path:
    pareto = singles[singles["measured_pareto_front"]].copy()
    path = FIG_DIR / "01_round1_lm_selected_measured_affinity_vs_stability.png"

    plt.figure(figsize=(11, 8))
    plt.scatter(singles["fitness"], singles["stability_deltaG"], s=42, alpha=0.38, color="#8a97a3",
                label="All fully labeled single mutants")
    plt.scatter(pareto["fitness"], pareto["stability_deltaG"], s=110, facecolors="none",
                edgecolors="#1f77b4", linewidth=1.8, label="Measured Pareto frontier")
    plt.scatter(selected["fitness"], selected["stability_deltaG"], s=90, color="#d62728",
                edgecolor="black", linewidth=0.6, label="Round-1 candidates selected by LM plausibility")
    add_wt_binding_annotation(singles["stability_deltaG"])
    for _, row in selected.head(8).iterrows():
        plt.annotate(row["mutation"], (row["fitness"], row["stability_deltaG"]), fontsize=8,
                     xytext=(4, 4), textcoords="offset points")
    plt.xlabel("Measured GB1 IgG-Fc binding fitness from database (higher is better)")
    plt.ylabel("Measured stability score from database (larger = more stable in source scale)")
    plt.title("Round 1 measured affinity/stability tradeoff: LM-plausible candidates only")
    plt.legend(frameon=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()
    return path


def plot_round1_validation(singles: pd.DataFrame, selected: pd.DataFrame) -> Path:
    pareto = singles[singles["measured_pareto_front"]].copy()
    path = FIG_DIR / "02_round1_measured_validation_no_lda.png"

    plt.figure(figsize=(11, 8))
    plt.scatter(singles["fitness"], singles["stability_deltaG"], s=42, alpha=0.35, color="#8a97a3",
                label="All fully labeled single mutants")
    plt.scatter(pareto["fitness"], pareto["stability_deltaG"], s=110, facecolors="none",
                edgecolors="#1f77b4", linewidth=1.8, label="Measured Pareto frontier")
    plt.scatter(selected["fitness"], selected["stability_deltaG"], s=90, color="#d62728",
                edgecolor="black", linewidth=0.6, label="LM-plausible Round-1 selections")
    add_wt_binding_annotation(singles["stability_deltaG"])
    plt.xlabel("Measured GB1 IgG-Fc binding fitness (higher is better)")
    plt.ylabel("Measured stability score from database (larger = more stable in source scale)")
    plt.title("Round 1 measured validation without LDA")
    plt.legend(frameon=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()
    return path


def plot_round2_validation(background: pd.DataFrame, highlighted: pd.DataFrame) -> Path:
    pareto = background[background["validation_pareto_front"]].copy()
    path = FIG_DIR / "04_round2_double_validation_no_lda.png"

    plt.figure(figsize=(11, 8))
    plt.scatter(
        background["measured_double_fitness"],
        background["double_stability_deltaG_additive"],
        s=14,
        alpha=0.12,
        color="#8a97a3",
        label="All double mutants with measured binding + additive stability",
    )
    plt.scatter(
        pareto["measured_double_fitness"],
        pareto["double_stability_deltaG_additive"],
        s=115,
        facecolors="none",
        edgecolors="#1f77b4",
        linewidth=1.8,
        label="Validation Pareto frontier",
    )
    plt.scatter(
        highlighted["measured_double_fitness"],
        highlighted["double_stability_deltaG_additive"],
        s=90,
        color="#d62728",
        edgecolor="black",
        linewidth=0.6,
        label="No-LDA doubles from LM-plausible singles",
    )
    add_wt_binding_annotation(background["double_stability_deltaG_additive"])
    for _, row in highlighted.head(5).iterrows():
        plt.annotate(row["double_mutation"], (row["measured_double_fitness"], row["double_stability_deltaG_additive"]),
                     fontsize=8, xytext=(4, 4), textcoords="offset points")
    plt.xlabel("Measured double-mutant GB1 binding fitness (higher is better)")
    plt.ylabel("Additive database stability score1 + score2 (larger = more stable in source scale)")
    plt.title("Round 2 double-mutant validation without LDA")
    plt.legend(frameon=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()
    return path


def add_wt_binding_annotation(y_values: pd.Series) -> None:
    plt.axvline(WT_BINDING_FITNESS, color="black", linestyle="--", linewidth=1.1, alpha=0.75)
    y_range = y_values.max() - y_values.min()
    y_bottom = y_values.min()
    plt.annotate(
        "WT binding threshold\nfitness = 1.0",
        xy=(WT_BINDING_FITNESS, y_bottom + 0.10 * y_range),
        xytext=(WT_BINDING_FITNESS + 0.18, y_bottom + 0.22 * y_range),
        arrowprops={"arrowstyle": "->", "color": "black", "lw": 1.0},
        fontsize=10,
        ha="left",
        va="bottom",
        bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "black", "alpha": 0.75},
    )


def write_readme(params: dict) -> None:
    readme = f"""# GB1 no-LDA co-optimization workflow

This workflow removes the LDA branch entirely.

## What the code does

1. Loads measured GB1 single-mutant binding and stability.
2. Treats local `ddg` as the source-reported stability score. This is not the lower-is-better folding reaction deltaG convention; larger values are treated as more stable in this dataset.
3. Scores each single mutant by LM plausibility proxy.
4. Selects the top LM-plausible single mutants.
5. Plots measured affinity/stability for selected singles.
6. Recombines selected singles into double mutants.
7. Looks up measured double-mutant binding.
8. Estimates double-mutant stability additively as `stability_score1 + stability_score2`.
9. Plots measured double-mutant binding versus additive stability.

## What was removed

No LDA models are trained.
No predicted LDA graph is generated.
There is no graph 3 in this package.

## Parameters

```json
{json.dumps(params, indent=2)}
```
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")


def main() -> None:
    ensure_clean_dirs()
    sns.set_theme(style="whitegrid", context="notebook")

    singles, checks, lm_score_source = load_and_score_singles()
    selected = select_round1_candidates(singles)
    doubles = pd.read_csv(DOUBLE_FITNESS_PATH)
    generated_doubles, highlighted_doubles = build_round2_candidates(selected, doubles)
    double_background = build_double_validation_background(singles, doubles)

    for src in [SINGLE_FITNESS_PATH, SINGLE_STABILITY_PATH, DOUBLE_FITNESS_PATH, STABILITY_SOURCE_PATH]:
        if src.exists():
            shutil.copy2(src, INPUT_DIR / src.name)

    plot_round1_clean(singles, selected)
    plot_round1_validation(singles, selected)
    plot_round2_validation(double_background, highlighted_doubles)

    params = {
        "run_created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "workflow": "no LDA",
        "wt_sequence_binding_table": WT_SEQUENCE_BINDING_TABLE,
        "wt_binding_fitness": WT_BINDING_FITNESS,
        "stability_metric": "local ddg column equals source deltaG/stability scale; larger values are treated as more stable in this dataset, not lower-is-better folding reaction deltaG",
        "lm_score_source": lm_score_source,
        "target_round1_candidates": TARGET_ROUND1_CANDIDATES,
        "actual_round1_candidates": int(len(selected)),
        "selected_round1_above_wt_binding": int((selected["fitness"] > WT_BINDING_FITNESS).sum()),
        "selected_round1_on_measured_pareto": int(selected["measured_pareto_front"].sum()),
        "selected_round1_mean_binding_fitness": float(selected["fitness"].mean()),
        "selected_round1_mean_stability_score": float(selected["stability_deltaG"].mean()),
        "single_mutants_with_binding_and_stability": int(len(singles)),
        "single_mutants_above_wt_binding": int((singles["fitness"] > WT_BINDING_FITNESS).sum()),
        "generated_round2_double_candidates": int(len(generated_doubles)),
        "generated_round2_doubles_with_measured_binding": int(generated_doubles["has_measured_double_fitness"].sum()),
        "double_validation_background_rows": int(len(double_background)),
        "double_validation_background_above_wt_binding": int(
            (double_background["measured_double_fitness"] > WT_BINDING_FITNESS).sum()
        ),
        "lda_removed": True,
        "graph_3_removed": True,
    }

    selected.to_csv(TABLE_DIR / "01_round1_lm_selected_candidates.csv", index=False)
    singles.sort_values("lm_score", ascending=False).to_csv(TABLE_DIR / "01_all_labeled_singles_scored.csv", index=False)
    checks.to_csv(TABLE_DIR / "stability_metric_crosscheck.csv", index=False)
    generated_doubles.to_csv(TABLE_DIR / "04_generated_round2_double_candidates.csv", index=False)
    highlighted_doubles.to_csv(TABLE_DIR / "04_highlighted_round2_doubles.csv", index=False)
    double_background.to_csv(TABLE_DIR / "04_double_validation_background.csv", index=False)
    pd.DataFrame([params]).to_csv(TABLE_DIR / "workflow_parameters.csv", index=False)
    (OUT / "workflow_parameters.json").write_text(json.dumps(params, indent=2), encoding="utf-8")

    manifest = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file():
            manifest.append({
                "relative_path": str(path.relative_to(OUT)),
                "size_bytes": path.stat().st_size,
            })
    pd.DataFrame(manifest).to_csv(OUT / "MANIFEST.csv", index=False)
    write_readme(params)

    print(json.dumps(params, indent=2))
    print(f"Results folder: {OUT}")


if __name__ == "__main__":
    main()
