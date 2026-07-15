# GB1 No-LDA Co-Optimization Demo

This repository contains a browser-based demo for a GB1 protein-variant co-optimization workflow.

The demo:

- scores single mutants using ESM/LM plausibility scores,
- selects Round-1 single-mutant candidates,
- validates selected candidates against measured GB1 binding and stability-source values,
- recombines selected singles into Round-2 double mutants,
- validates generated doubles against measured double-mutant binding data.

No LDA model is used in this current workflow, and there is intentionally no Graph 3.

## Website

Open `index.html` in a browser, or use the GitHub Pages link if Pages is enabled for this repository.

## Input files

Use the files in `data/`:

| Website input | File |
|---|---|
| Single-mutant binding CSV | `data/gb1_fitness_singles.csv` |
| Single-mutant stability CSV | `data/gb1_stability_deltaG.csv` |
| Double-mutant binding CSV | `data/gb1_fitness_doubles.csv` |
| Optional ESM/LM scores CSV | `data/gb1_lm_scores.csv` |

## Important interpretation notes

- `lm_score` is a z-scored ESM WT-marginal log-likelihood-ratio. It is a sequence-plausibility score, not measured affinity or stability.
- Graph x/y positions use measured database values.
- The plotted stability value is the source-reported GB1 stability scale used in this project. It should not be relabeled as generic folding deltaG where lower/more negative is always favorable.
- Double-mutant stability is estimated additively as `stability_score_1 + stability_score_2`.

## Local workflow files

The core workflow used to generate the static outputs is included:

- `run_gb1_no_lda_workflow.py`
- `score_gb1_esm.py`

