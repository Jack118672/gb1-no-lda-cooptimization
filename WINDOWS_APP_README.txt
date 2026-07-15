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
- Automatically loads a bundled real GB1 WT PDB structure, so no extra structure files are required.
- Lets you inspect selected red or blue graph points, map mutation sites onto the real GB1 structure trace, and view an estimated local mutant shape effect.
- Optionally lets you upload WT and mutant PDB files from AlphaFold, Rosetta, molecular dynamics, docking, ColabFold, or experimental structures if you want to override the bundled model.

Important:
The 3D Mutation Model uses a bundled real GB1 WT PDB structure by default. The mutant view is generated automatically by locally deforming that WT structure around selected mutation sites. This is useful for visualization, but it is still not molecular dynamics, docking, AlphaFold, Rosetta, or a true new folding prediction.
For true predicted mutant structures, run the prediction tool externally, export WT and mutant PDB files, then upload those PDB files in the 3D Model tab.
