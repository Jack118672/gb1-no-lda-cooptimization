# GB1 Multi-Mutant Co-Optimization App

Browser app for a no-LDA GB1 protein-variant co-optimization workflow.

Live site:
https://jack118672.github.io/gb1-no-lda-cooptimization/

Windows app download:
https://jack118672.github.io/gb1-no-lda-cooptimization/downloads/GB1_Multi_Mutant_App_Windows.zip

## What The App Does

- Loads packaged GB1 binding, stability, double-mutant, and ESM/LM score CSV files.
- Scores single mutants using ESM/LM plausibility.
- Selects Round-1 single-mutant candidates.
- Validates selected candidates against measured GB1 binding and source-reported stability values.
- Recombines selected singles into double mutants and checks measured double-mutant binding.
- Extends beyond doubles into capped higher-order multi-mutant recommendations, now up to 8 mutations.
- Provides an interactive 3D Mutation Model that automatically loads a bundled real GB1 WT PDB structure, renders it as a protein-cartoon ribbon with helix, beta-sheet, and loop regions, maps selected mutation sites onto that structure, and shows an estimated local mutant shape effect.

No LDA model is used in this current workflow.

## Use In A Browser

Open the live site:
https://jack118672.github.io/gb1-no-lda-cooptimization/

Then:

1. Open **1. Inputs / Run**.
2. Click **Run**.
3. Open **2. Outputs / Results** to view graphs, tables, and summaries.
4. Click red or blue graph points to open **3. 3D Model**.
5. The **3D Model** tab automatically uses a bundled real GB1 WT PDB structure and shows the before/after mutation comparison.

## Download As A Windows App

Download:
https://jack118672.github.io/gb1-no-lda-cooptimization/downloads/GB1_Multi_Mutant_App_Windows.zip

Then:

1. Extract the zip.
2. Double-click **GB1 Multi-Mutant App.bat**.
3. A browser window opens automatically.
4. Keep the launcher window open while using the app.
5. Close the launcher window when finished.

The Windows launcher uses PowerShell to serve the app locally, so the browser can load the packaged CSV files.

## Input Files

The app uses the files in `data/`:

| App input | File |
|---|---|
| Single-mutant binding CSV | `data/gb1_fitness_singles.csv` |
| Single-mutant stability CSV | `data/gb1_stability_deltaG.csv` |
| Double-mutant binding CSV | `data/gb1_fitness_doubles.csv` |
| ESM/LM scores CSV | `data/gb1_lm_scores.csv` |

## Important Interpretation Notes

- `lm_score` is a z-scored ESM WT-marginal log-likelihood-ratio. It is a sequence-plausibility score, not measured affinity or stability.
- Graph x/y positions use measured database values when available.
- The plotted stability value is the source-reported GB1 stability scale used in this project. It should not be relabeled as generic folding deltaG where lower/more negative is always favorable.
- Double-mutant stability is estimated additively as `stability_score_1 + stability_score_2`.
- Higher-order multi-mutants beyond measured doubles are ranked recommendations, not experimentally measured binding claims in the local dataset.
- The 3D Mutation Model uses a bundled real GB1 WT PDB structure by default. It highlights mutation positions on that structure, compares WT vs mutated sequence, and automatically generates a local mutant shape effect around changed residues. It does not run molecular dynamics, docking, AlphaFold, Rosetta, or a true new folding prediction.
