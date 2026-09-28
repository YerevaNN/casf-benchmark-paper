# Manuscript figures

Four figures are included through `../figures.tex` and `../results.tex`:

1. `core-recovery-budget`: core recovery at 0.75 Å and retained conformer counts
   at both candidate targets, for all ten pipelines and a stored reference.
2. `diversity-recovery`: all ten pipelines plus ChEMBL3D-PB, core recovery versus
   cluster count at the separate 1.0 Å clustering radius.
3. `size-flexibility`: core-only heatmaps with recovered/total counts for every
   method in the main panel; groups with fewer than five ligands are marked.
4. `multireference-coverage`: exploratory multi-reference recall/precision and
   paired molecules; coincident points show multiplicity.

Each has a vector PDF, editable SVG, and 450-dpi PNG at 7 inches wide.
`main-figures.pdf` is regenerated with the same four figures. Sources and
provenance are archived in `../figure-data/`.

Run `python jcim-overleaf/build_figures.py` from the paper repository root in
an environment containing NumPy, pandas, and Matplotlib. Figure rendering does
not compile or update manuscript previews.

`energy-window-historical` is unchanged historical evidence, excluded from the
paper and gallery. Rendering it requires the explicit `--historical-energy`
option. The removed multi-threshold figure remains available in Git history.
