GB1 Multi-Mutant App

How to launch:
1. Extract this folder from the zip file.
2. Double-click "GB1 Multi-Mutant App.bat".
3. A browser window opens automatically.
4. Keep the black launcher window open while using the app.
5. Close the launcher window when you are done.

What the app does:
- Loads the packaged GB1 binding, stability, double-mutant, and ESM/LM score CSV files.
- Runs the no-LDA co-optimization workflow.
- Shows Inputs / Run, Outputs / Results, and a 3D Mutation Model.
- Lets you inspect selected red or blue graph points, map mutation sites onto a simplified GB1 fold, and view an estimated local shape perturbation.
- Lets you upload WT and mutant PDB files from AlphaFold, Rosetta, molecular dynamics, docking, ColabFold, or experimental structures and view the actual coordinate traces side by side.

Important:
The 3D Mutation Model is an interpretive structure map. It highlights mutation positions on a simplified GB1 fold, compares WT vs mutated sequence, and bends the mutant-side view near the changed residues using amino-acid chemistry differences. It does not run molecular dynamics, docking, AlphaFold, Rosetta, or a true new folding prediction.
For real predicted structures, run the prediction tool externally, export WT and mutant PDB files, then upload those PDB files in the 3D Model tab.
