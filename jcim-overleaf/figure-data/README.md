# Figure data from the public dashboard

Source: [dashboard-data-stereo-identity-v2](https://github.com/YerevaNN/casf-benchmark/releases/tag/dashboard-data-stereo-identity-v2),
published 18 September 2026 and still GitHub's latest release when checked on
28 September 2026. The public dashboard mirror's default loader fetches these
assets from YerevaNN/casf-benchmark. Both local database SHA-256 hashes exactly
match GitHub's asset digests; see `public-release.json` and `provenance.json`.
This identifies the public data release, not the cache state of an individual
Streamlit browser session.

## Files and denominators

- `casf_summary.csv`: plotted recovery thresholds, cluster radii, candidate tiers,
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

Main figures use all 94 core or 1,236 ref entries for recovery, including
missing outputs as misses. Clustering means are conditional on defined values.
ChEMBL3D-PB is a stored ensemble, not a newly generated 1,000-sample baseline.
Sampling in Figure 1B is a comparison of candidate-tier endpoints; the graph
is not a random fixed-valid-K curve.

Energy results in the release are historical and uncorrected relative to the
main database. The diagnostic identifies that status visibly and is excluded
from the manuscript. No new force-field energies or conformers were generated.

## Reproduction

Run `python manuscript/build_figures.py` in the project environment to render
from the archived CSVs. Add `--refresh-data` to re-export from local databases;
the script requires their hashes to match the pinned public release metadata.
The renderer uses pandas, NumPy, and Matplotlib. The 7-inch vector PDFs embed
TrueType fonts; SVGs retain editable text, and color PNGs are exported at 450 dpi.
The method palette is consistent across figures, with symbols/line styles as
additional identifiers. Figure panels intentionally omit additional Qwen
checkpoints and minimized classical variants; Table 1 retains those classical
comparisons, while Figure 3 uses the four contrasting methods in Table 2.
