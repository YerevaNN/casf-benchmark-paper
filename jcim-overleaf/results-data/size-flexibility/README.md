# Figure 5: size and flexibility across three datasets

Run `python jcim-overleaf/results-data/size-flexibility/build_figure.py` from the
paper repository root (requires NumPy, pandas, and matplotlib). It rebuilds
`figures/figure-4-size-flexibility.{pdf,svg,png}` and the derived CSVs here.
The filename is retained for stable manuscript links; this is main Figure 5.

Rows are PLINDER-22, CASF core94, and CASF ref1236. Columns are heavy atoms and
rotatable bonds. Method order follows Table 1. The original Qwen 1.7B FSQ
step47023 is retained; these are not retrained results.

- `casf-size-source.csv`: unchanged copy of benchmark
  `docs/publication_tables_2026_09_28/size_strata.csv`.
- `casf-flexibility-source.csv`: unchanged copy of benchmark
  `docs/publication_tables_2026_09_24/flexibility.csv`.
- `source-provenance.json`: original paths and hashes for those archives.
- `figure5-strata.csv`: plotted values, metrics, group sizes, and CASF hit counts.
- `plinder-strata-per-molecule.csv`: PLINDER identities, descriptors, groups, and COV-R.
- `provenance.json`: hashes of all current inputs and validation summary.

CASF uses Hit@0.75 from PB-filtered fixed-target pools; the stored baseline uses
its original pool. The builder reproduces every one of the 88 previous core
cells in `../strata.csv` and verifies group totals of 94/1236 for each method.

PLINDER joins `../../plinder22-data/per-molecule.csv` to the cleaned inventory in
`../../method-data/plinder22-molecules.csv` by exact SMILES. Each method has the
same 22 identities and finite coverage. Stratum COV-R is an equal-weight mean
over molecules, not a pooled reference fraction or binary molecule hit rate.
The supplied pools are not PB-filtered. Reference matching uses RMSD <=0.75 Å.

Heavy-atom bins are <20, 20–29, 30–39, >=40; rotatable-bond bins are 0–3, 4–6,
7–8, >=9. Asterisks/dashed outlines mark n<5; empty groups and unavailable
methods are gray dashes. PLINDER has no >=40-heavy-atom molecules. Pool
conditions and outcome definitions differ across datasets; the shared color
scale does not establish direct performance comparability. No conformers,
RMSD alignments, or validity checks are recalculated by this renderer.
