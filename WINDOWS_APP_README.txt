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
- Generates capped higher-order multi-mutant recommendations up to 8 mutations.
- Shows Inputs / Run, Outputs / Results, and a 3D Mutation Model.
- Automatically loads a bundled real GB1 WT PDB structure, so no extra structure files are required.
- Draws the 3D Model as a protein-cartoon ribbon with helix, beta-sheet, and loop regions instead of a dot cloud.
- Lets you inspect selected red or blue graph points, map mutation sites onto the real GB1 cartoon structure, and view an estimated local mutant shape effect.

Important:
The 3D Mutation Model uses a bundled real GB1 WT PDB structure by default. The mutant view is generated automatically by locally deforming that WT structure around selected mutation sites. This is useful for visualization, but it is still not molecular dynamics, docking, AlphaFold, Rosetta, or a true new folding prediction.
