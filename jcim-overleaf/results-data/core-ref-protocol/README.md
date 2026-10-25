# Completed core/ref energy and sampling figures

This archive supplies Figures 3, 4, and B1 and Supporting Table S9, integrated
on 26 October 2026 from the approved completed previews. Figure colors, method
orders, paired dataset layout, scales, and values are preserved; review-only
preview subtitles are removed. Main Tables 1–4 are unchanged.

## Reproduction

From the paper repository, with NumPy, pandas, and Matplotlib installed:

```bash
python jcim-overleaf/results-data/core-ref-protocol/build_figures.py
python jcim-overleaf/results-data/core-ref-protocol/build_energy_table.py
```

These commands render only the figure artwork and S9 source, not the manuscript.
The three figures retain their existing PDF/PNG/SVG filenames under `figures/`.
The renderer reads only files in the paper repository. No database access,
force-field calculation, new sampling, or network access is needed to render.
`provenance.json` records the manuscript-relative input paths and SHA-256 hashes.
`approved-preview-provenance.json` retains the original preview input hashes,
including the independently checked public-dashboard database hash.

## Inputs and calculation lineage

- Core energies: `../../energy-data/summary.csv` and `protocol.json`.
- Core curves: `../../sampling-data/`, including display order and palette.
- Two-target endpoints: `../three-dataset-tables/table3.csv` for both datasets.
- Ref energy summaries and per-entry audit: `ref/summary.csv` and
  `ref/per_molecule.csv`; the latter retains the failed-energy count.
- Ref curves, 80% frontier, and stored-ensemble comparison: `ref/recovery_curves.csv`,
  `ref/recovery80_frontier.csv`, and `ref/reference.json`.
- Ref full-pool endpoints independently checked against the published database:
  `ref/full-pool-endpoints.csv`.
- Calculation settings and original source hashes: `ref/provenance.json`.

The completed calculation code and results were committed in benchmark commit
[`5090d110671731fefbe31c81fb65258c178e6304`](https://github.com/YerevaNN/casf-benchmark/commit/5090d110671731fefbe31c81fb65258c178e6304).
The calculation is implemented in `docs/ref_protocol_analysis/calculate.py` and
`src/casf_benchmark/analysis/common_protocol.py` in that repository. The
[public data release](https://github.com/YerevaNN/casf-benchmark/releases/tag/dashboard-data-ref-protocol-v1)
includes `core-ref-protocol-results.tar.gz` with the full per-conformer results,
candidate distances, audit records, and core inputs. Those large files are not
duplicated in the paper repository. The dashboard SQLite SHA-256 used for the
independent endpoint checks is
`626ac32e96a302bfb0e8c6c9c49a9f529c023366ac193b5805b5892cd84be804`.

All panels retain Qwen 1.7B FSQ step47023. These are not retrained-model results.

## Protocol and completeness

Both datasets use RDKit 2025.09.6 and MMFF94s, removal and rebuilding of hydrogens
with assigned coordinates, no coordinate optimization, nonbonded threshold 100,
and ignored interfragment interactions. Population energy SD uses `ddof=0`.
Heavy-atom coordinates are preserved to 1e-10 Å. Energy figures pair cluster means
and median energy SD on the same measurable entries for each method.

Saved PoseBusters membership is preserved. The ref stored comparison was
reconstructed from the original stereochemically matched coordinates; archived
energies verified conformer identity and order. Six filtered entries used uniquely
matched ordered energy records to recover the original retained subset. Replaying
the stochastic PoseBusters energy test can change that subset and was not used
for the final energy pools. All energies shown here are newly evaluated under
the common protocol, rather than copied from those identity-check records.

Core has 77,986 finite energies and three empty method–entry outputs. Ref has
835,972 finite energies out of 835,973 retained conformers and 50 empty outputs.
The one undefined ref energy is in FlowR entry `2qm9_2qm9_conf0`; it is omitted
from its ensemble's energy moments and counted explicitly. Finite outliers,
including the large ref LoQI mean SD, remain in the summaries. Scientific notation
in S9 shortens the printed value without changing the calculation.

Sampling uses exact uniform subsets of the saved pre-PoseBusters candidate pools,
with rejected candidates consuming budget and empty retained pools contributing
zero recovery. Budgets are capped at available candidate counts. Each dataset
uses all its entries with equal weight: 94 core or 1,236 ref. There are 63 distinct
budgets from 10 to 1,000; the frontier inverts expected recovery over observed RMSD
values. Confidence intervals from 2,000 paired bootstrap resamples (seed 20261006)
are archived but not drawn. Ref contains 12,351,409 candidates and 11,597,494
retained conformers with finite reference RMSDs. Its full-pool endpoints agree
with the published database; the core protocol reproduces archived curves and
intervals within 1e-14 in the upstream parity audit.

The completed figure values are checked against the approved previews. S9 core
values retain their prior rounding; ref rows are aggregated from the archived
per-entry results. Source citation/cross-reference checks accompany the edit.
No local manuscript compilation is performed; Overleaf is the review target.
