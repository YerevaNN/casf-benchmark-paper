# Figure 5: size and flexibility in CASF core and ref

Run `python jcim-overleaf/results-data/size-flexibility/build_figure.py` from the
paper repository root (requires NumPy, pandas, and matplotlib). It rebuilds
`figures/figure-4-size-flexibility.{pdf,svg,png}` and `figure5-strata.csv` here.
The filename is retained for stable manuscript links; this is main Figure 5.

Rows are CASF core94 and CASF ref1044. Columns are heavy atoms and
rotatable bonds. Method order follows Table 1. The original Qwen 1.7B FSQ
step47023 is retained; these are not retrained results.

- `casf-size-source.csv`: unchanged core rows from benchmark
  `docs/publication_tables_2026_09_28/size_strata.csv`; ref rows are recomputed from
  the selected per-entry scores archived in `../ref-panel/`.
- `casf-flexibility-source.csv`: unchanged core rows from benchmark
  `docs/publication_tables_2026_09_24/flexibility.csv`; ref rows use the same selection.
- `source-provenance.json`: historical original paths and hashes. Current ref
  selection and source hashes are in `../ref-panel/provenance.json`.
- `figure5-strata.csv`: plotted values, metrics, group sizes, and CASF hit counts.
- `plinder-strata-per-molecule.csv`: historical PLINDER breakdown, retained as an archive and excluded from the current figure.
- `provenance.json`: hashes of all current inputs and validation summary.

CASF uses Hit@0.75 from PB-filtered fixed-target pools; the stored baseline uses
its original pool. The builder reproduces every one of the 88 previous core
cells in `../strata.csv` and verifies group totals of 94/1044 for each method.

Heavy-atom bins are <20, 20–29, 30–39, >=40; rotatable-bond bins are 0–3, 4–6,
7–8, >=9. Asterisks/dashed outlines mark n<5; empty groups are gray dashes.
The original palette and method order are preserved. No conformers, RMSD
alignments, or validity checks are recalculated by this renderer.
