# Additional Qwen evaluation

Run: `qwen_0p6b_casf_dedup_binned_pm12_step64000` (0.6B, binned coordinates).
Source: the read-only `plinder_per_molecule` and `plinder_summary` tables in
`/mnt/weka/vtarasov/outputs/casf_benchmark_runs/dashboard_local/extended_casf_analysis.sqlite`.
The 22 molecule identities and 1,067 reference counts match the paper panel.
All original method scores remain unchanged.

`per-molecule.csv` and `summary.csv` preserve the source records.
`table-records.csv` changes only the 0.75 coverage column names to the existing
Table 4 convention. The table renderer reads this supplement alongside the
original ten-method records. `nearest-distances.npz` and `entries.json` archive
both directions of nearest distances from the existing RMSD cache; all minima
are finite. Four nonfinite pairwise comparisons have finite alternatives and
do not exclude any reference or generated conformer. No new scoring was run.

Run `python jcim-overleaf/results-data/build_druglike_table.py` to render Tables
4 and S4. Run `python jcim-overleaf/plinder22-data/qwen-dedup/verify.py` for the
supplemental scores and figure checks, and the parent `verify.py` for the original
panel, source references, and table reproducibility.

The run label identifies the training recipe; the actual exclusion manifest is
not supplied here. Do not infer PLINDER training exclusion from the label or
interpret the comparison with the 1.7B FSQ model as a controlled deduplication ablation.
