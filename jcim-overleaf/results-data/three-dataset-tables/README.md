# Current three-dataset table sources

## Current ref data, 30 October 2026

All ref rows now use the 1,044-crystal panel from the published `ref_dedup`
cohort. The paper calls it ref. Core and PLINDER rows are unchanged. Current
ref inputs, the one-crystal-per-identity selection manifest, delta audit, and
numerical checks are in `../ref-panel/`. Rebuild Tables 1–3 with `build_tables.py`
and Figure 2 with `plot_casf.py`. The older source descriptions below record
previous exports and are superseded for ref by this update.


## Refined wide layout, 22 October 2026

The author selected option A: shorter headers, full-width tabular spacing,
light row shading, and concise notes. The same renderer preserves all source
values, method orders, and generated-only best/second-best emphasis. Tables
1/2 remain PLINDER/core/ref; Table 3 remains CASF-only with paired targets.

PLINDER rows in Tables 1 and 2 were replaced on 21 October 2026 using the
validated cleaned release in `../../plinder22-data/`. CASF rows and Table 3
are unchanged. `plinder-per-molecule.csv` contains the refreshed source records.
The original preview provenance files and cache/table comparison below describe
the superseded PLINDER-23 snapshot and remain historical CASF/source lineage.
Current PLINDER hashes and validation are in `../../plinder22-data/provenance.json`.
Rebuild with `build_tables.py`; `plot_casf.py` renders the two-panel CASF scatter.

## Original preview archive (historical PLINDER values)

# Tables 1–3: approved wide layouts

Integrated 19 October 2026 from the author's reviewed previews in the analysis
repository. Table 1 and Table 2 use PLINDER/core/ref column groups. Table 3 has
core/ref groups and paired candidate-target rows; PLINDER was explicitly removed.

Run `python jcim-overleaf/results-data/three-dataset-tables/build_tables.py` from
the paper repository root to rebuild `diversity.tex`, `count-recovery.tex`, and
`recovery.tex`. The renderer uses the standard library and these CSVs only. It
preserves existing row orders and ranks distinct unrounded values within each
dataset (and budget for Table 3), excluding the stored ChEMBL3D baseline.

- `table1.csv`: clustering means, missing-value status, and defined denominators.
- `table2.csv`: reference coverage/Hit percentages, matching distances, denominators.
- `table3.csv`: CASF budget comparisons, denominators, and mean retained counts.
- `plinder-per-molecule.csv`: reference counts and COV-R/MAT-R used in Table 2.
- `cache-table-comparison.csv`: Table 2 cache values versus earlier Table 4 records.
- `table*-provenance.json`: unchanged provenance snapshots from the approved
  previews, including original input paths/hashes and aggregation rules. Their
  “manuscript unchanged” fields describe the preview stage, before integration.

CASF sources use the database hash recorded in `../sampling-budget-provenance.json`.
Read-only extraction selects `mol_id`, `casf_best_rmsd`, `greedy_clusters_0p5`,
`greedy_clusters_1p0`, `largest_cluster_fraction_1p0`, and `post_pb_confs` from
`per_ligand_long`, restricted by `ligand_set` and the selected method IDs.
Hit uses all entries; other metrics use defined entries. Core means were verified
against the existing manuscript exports. Reference records are complex entries,
not 1,236 distinct parent compounds.

PLINDER clustering comes from `../druglike_selected.csv` and uses up-to-200
supplied-conformer subsets. Missing largest-cluster fractions stay pending.
Table 2 uses cached reference minima for all eight methods, strict RMSD cutoffs,
and equal molecule weights. All reference minima are finite. The original six
learned-method nearest distances are in `../../plinder-threshold-data/`; the two
classical cache hashes are recorded in Table 2 provenance. No alignment,
clustering, generation or force-field calculation was rerun.

All values remain PLINDER-23 and the original Qwen 1.7B FSQ step47023 checkpoint.
They do not represent the ongoing retraining or the planned PLINDER-22 panel.
PLINDER supplied-pool filtering differs from CASF. Small cache/Table 4 differences
remain unresolved and Table 4 has not been overwritten.
