# Figure data from the public dashboard

Source: [dashboard-data-stereo-identity-v2](https://github.com/YerevaNN/casf-benchmark/releases/tag/dashboard-data-stereo-identity-v2),
published 18 September 2026 and still GitHub's latest release when checked on
28 September 2026. The public dashboard mirror's default loader fetches these
assets from YerevaNN/casf-benchmark. Both local database SHA-256 hashes exactly
match GitHub's asset digests; see `public-release.json` and `provenance.json`.
This identifies the public data release, not the cache state of an individual
Streamlit browser session.

## Files and denominators

- `core-recovery-counts.csv` and `core-recovery-audit.csv`: independent recount
  of the five opening Results comparisons, with all 94 ligand records per
  method and the minimum RMSD used to classify recovery.

- `casf_summary.csv`: archived recovery thresholds, cluster radii, candidate tiers,
  and retained counts, calculated from the public database's per-ligand rows.
- `casf_selected_entries.csv`: selected fixed-tier source entries, including
  common reference descriptors; these reproduce the size/flexibility groups.
- `casf_strata.csv`: group sizes and full-group recovery for those descriptors.
- `dashboard_displayed_aggregates.csv`: public summary rows for auditing. Some
  dashboard percentages omit undefined outcomes. For example, core MCF and
  NExT-Mol show 81/93 = 87.1% among measurable entries; the manuscript uses
  81/94 = 86.2% over all mapped entries. This is an explicit denominator change,
  not a different set of generated samples.
- `drug_molecule_metrics.csv`: 23 molecule rows for each of the six selected
  drug methods. These are supplied-pool metrics without a common PB filter.
- `drug_coverage_intervals.csv`: means and marginal 95% molecule-bootstrap
  intervals, 10,000 resamples, NumPy seed 20260924. They are not the paired
  method-difference intervals quoted in the Results.
- `historical_energy_summary.csv`, `historical_energy_plotted.csv`: unchanged
  14 September sidecar and the subset plotted in the separate diagnostic.

Main Qwen identity is explicitly `qwen_1p7b_fsq_bigdata_step47023`. The public
release also contains legacy rows with names such as
`qwen_1p7b_fsq_bigdata_pretrain`; they are not interchangeable.

Main CASF figures use all 94 core entries for recovery; supporting tables use
1,236 ref entries. Both include
missing outputs as misses. Clustering means are conditional on defined values.
ChEMBL3D-PB is a stored ensemble, not a newly generated 1,000-sample baseline.
Figure 1 compares recovery at the two candidate-tier endpoints; it
is not a random fixed-valid-K curve.

Energy results in the release are historical and uncorrected relative to the
main database. The diagnostic identifies that status visibly and is excluded
from the manuscript. No new force-field energies or conformers were generated.

## Reproduction

Run `python jcim-overleaf/build_figures.py` in the project environment to render
from the archived CSVs. Add `--refresh-data` to re-export from local databases;
the script requires their hashes to match the pinned public release metadata.
The renderer uses pandas, NumPy, and Matplotlib. The 7-inch vector PDFs embed
TrueType fonts; SVGs retain editable text, and color PNGs are exported at 450 dpi.
The method palette is consistent across figures, with symbols/line styles as
additional identifiers. Figures 1--3 include all ten main pipelines, including the minimized classical
variants, with one representative Qwen checkpoint. Ref is reported in SI.

## Recovery tables at both sampling targets

`recovery-table-tiers.csv` archives the 42 selected cohort/method/tier rows
used in Table 1 and supporting Tables S1 and S2. It is a subset of
`docs/publication_tables_2026_09_24/all_casf_summary.csv` in the benchmark
repository; hit rates are converted from fractions to percentages. Recovery
and retained counts were checked against `casf_summary.csv` for every
overlapping row. The subset also includes the minimized classical methods.
RMSD averages use entries with defined values; recovery and mean retained
counts use all mapped entries, including missing outputs. The stored
ChEMBL3D-PB ensemble is a single reference, not a separately generated tier.

## Core-only figures at 0.75 Å

- `core-comparison.csv`: 11 fixed/reference rows for all main methods, with
  mean retained count, mean cluster count at 1.0 Å, and recovery at 0.75 Å.
- `core-strata.csv`: 88 core rows (11 methods × 2 descriptors × 4 groups),
  with exact recovered/total counts. All overlapping values agree with the
  earlier `casf_strata.csv`. Integer counts reconstructed from the flexibility
  export were checked against its percentages.
- `core-figure-sources.json`: original corrected export paths, SHA-256 digests,
  checkpoint identity, and transformations for these new archived subsets.

The single recovery cutoff does not change the distinct clustering definition.
The added Table 2 mean supplied counts come from `drug_molecule_metrics.csv`.
Its MAT columns were moved to Table S4 without changing their values. Figure 4
preserves the original bootstrap draw order despite reordered display rows.
