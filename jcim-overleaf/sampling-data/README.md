# Combined sampling-budget figure

The displayed figure combines two panels for Qwen, LoQI, FlowR, and NExT-Mol.
Panel A shows only Hit@0.75; panel B shows the RMSD threshold for 80% expected
recovery. The archived analysis retains all ten fixed-target pipelines and
thresholds on the same core94 panel as main Table 3. Qwen is the initial 1.7B FSQ step47023 model;
benchmark-excluded retraining is still pending. These are new calculations on
the saved pools, not new generation runs or reuse of the older valid-only K curves.

## Candidate populations and RMSDs

`extract.py` reads the exact fixed-target SDFs identified by the pinned dashboard
snapshot and recomputes heavy-atom symmetry-aware aligned RMSDs with the existing
metric. All 884,399 retained conformers have finite RMSDs. Every per-entry minimum
and median agrees with the database. No coordinates are optimized. All RMSD
thresholds use the inclusive <= convention of the CASF benchmark.

The 940 method–entry pools contained 939,646 candidates before PoseBusters,
of which 55,247 were rejected. Rejections remain non-hits in each sampling
population. Two pools have no retained structures and remain in the denominator
of 94. Seventeen pools have fewer than 1,000 candidates (minimum 725); requested
budgets are capped at the available pre-PoseBusters count. Shortfalls are not
filled with invented conformers. Budget does not measure raw generator attempts,
upstream preparation failures, wall time, or GPU cost.

The compressed `candidate_rmsds.npz` contains sorted retained RMSDs padded with
infinity, the actual candidate counts, method/entry identifiers, and stored
ChEMBL3D per-entry minima. Padding beyond the actual count is never sampled.
Order is irrelevant for uniform subset sampling; these arrays cannot reproduce
first-K generation-order curves. `pool_audit.csv` records SDF paths/hashes,
counts, calculation fingerprints and failures. The source database, scripts,
RDKit version, and reference-file hashes are in the provenance records.

## Exact recovery calculation

For a pool of N candidates, h retained structures with RMSD <= t, and
k = min(requested budget, N), the probability of at least one hit in a uniformly
sampled subset without replacement is

    p(k,t) = 1 - choose(N-h,k) / choose(N,k).

`analyze.py` evaluates this hypergeometric probability directly rather than
introducing Monte Carlo noise from repeated subset draws. Expected recovery is
the equal-weight mean of the 94 probabilities. All 1,000-candidate Hit@0.5 and
Hit@0.75 endpoints reproduce the current manuscript; Hit@1.0 is calculated from
the same RMSDs. Budgets span 10–1,000 on an integer logarithmic grid including
10, 25, 50, 100, 250, 500 and 1,000. Lines join evaluated budgets.

The archive retains pointwise 95% percentile intervals from 2,000 bootstrap
resamples of the 94 entries, with seed 20261006. The identical entry weights are
used across methods, thresholds, and budgets. These intervals describe molecule
composition uncertainty conditional on the saved pools, not variation across
independent model runs or new generation pools. They are not simultaneous bands,
and overlap alone does not test a paired method difference. The approved figure
does not display these intervals. Both budget axes are linear, panel A uses a
50–95% recovery range, and endpoint values are labeled directly.

Panel B inverts expected recovery: at each budget it finds the
smallest observed RMSD threshold with mean recovery at least 80%. Binary search
uses the actual observed RMSDs, not a discretized threshold grid. This is the
threshold of the mean recovery function, not a mean of thresholds computed for
random draws. The 80% target is a presentation choice; this panel is descriptive
and has no uncertainty bands. At the full pool it agrees with the 76th ordered
per-entry minimum (at least 80% of 94 entries).

ChEMBL3D is always its full stored ensemble. Its dashed horizontal lines are
fixed reference levels, not a budget-matched or subsampled ChEMBL3D curve.
No curve is extrapolated beyond the saved pools.

## Outputs and reproduction

- `figure-budget-combined.*`: the current four-method, two-panel display.
- `display.json`: selected methods and the 0.75 Å recovery cutoff for panel A.
- The original separate figures remain historical assets; they are no longer included
  in the manuscript. Their renderer is available in the earlier repository history.
- `recovery_curves.csv`, `recovery80_frontier.csv`: exact plotted values.
- `entry_probabilities.npz`: per-entry probabilities for all plotted thresholds.
- `endpoints.csv`, `reference.json`: endpoint and stored-reference checks.
- `config.json`, `palette.json`: calculation settings and shared method colors.

Extraction assumes the original analysis checkout and Weka input paths:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python docs/sampling_analysis/extract.py
OPENBLAS_NUM_THREADS=1 python docs/sampling_analysis/analyze.py
python docs/sampling_analysis/plot.py
```

Once extracted, `analyze.py` and `plot.py` can run directly beside the archived
NPZ/config/palette files in the manuscript repository. Source SDFs and the main
benchmark database are never modified. `cache/` is a resumable execution aid,
not part of the manuscript archive. The combined figure uses the established method colors. Only the stored
ChEMBL3D reference levels are dashed; no minimized variants are displayed.

Validation checked all current-table endpoints, within-pool minima/medians,
monotonic recovery and thresholds, the frontier's crossing on both sides,
and the probability formula against direct combinations and random subsets.
