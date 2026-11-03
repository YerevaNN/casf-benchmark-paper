# Primary Qwen FSQ update

The primary Qwen model is `qwen_1p7b_v2_coordtoken_step41054`, displayed as
**Qwen 1.7B FSQ**. It replaces `qwen_1p7b_fsq_bigdata_step47023` throughout the
active paper, including the optional review appendix. The separate 0.6B
CASF-dedup model is preserved. Other methods and reference panels are unchanged.

## Sources and scope

- CASF results: the published `casf_analysis_dashboard.sqlite`, restricted to
  this exact family, core94 and the frozen `ref_dedup` selection of 1,044 crystals.
  The manuscript calls the latter **ref**. Missing distances count as recovery
  misses; continuous summaries use measurable entries. The new Qwen count target
  has 1,043 measurable ref entries, while its fixed target has all 1,044.
- PLINDER: the run's published per-molecule and summary files, checked against
  its cached RMSD matrices. The same 22 molecules and 1,067 references are used.
  Full supplied ensembles determine COV/MAT and clustering; PB is separate.
- Source hashes, inference metadata, coordinate-decoder identity, and the model
  export identity are archived here. The author identifies this run as the latest
  FSQ training. An executed training exclusion manifest remains outstanding;
  the run label is not evidence of a leakage-free evaluation.

## Energy and sampling

The publication did not include paper-protocol energy or sampling data for this
checkpoint. `calculate.py` therefore measures the saved PB-passing SDFs without
rerunning generation, filtering, or optimization. Explicit hydrogens are rebuilt
and MMFF94s single-point energies are evaluated with the original paper settings.
Heavy-atom coordinates are checked for zero displacement. Finite outliers remain.
The within-ensemble SD uses ddof=0. Current totals are 78,058 finite core energies
and 734,951 finite ref energies out of 734,952 retained conformers. The one
undefined ref energy remains the existing FlowR case.

Fixed-pool symmetry-aware aligned RMSDs are checked against the published
per-entry minima and medians. Exact recovery probabilities use uniform subsets
without replacement, with rejected candidates retained in the candidate budget.
The existing 63 budgets, thresholds, paired bootstrap seed and 80% frontier
inversion are unchanged. `core/` and `ref/` contain the new distance arrays,
per-conformer energies, aggregates, and curves; `calculation-provenance.json`
records every source SDF hash. External caches only accelerate reproduction.
SDFs are read sequentially, avoiding an RDKit list-length pre-scan issue without
changing any molecule or dropping any record.

## Reproduction

From `jcim-overleaf/`, with the benchmark environment available:

```bash
python results-data/qwen-v2/update.py
PYTHONPATH=/mnt/weka/mbedrosian/code/casf-benchmark/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python results-data/qwen-v2/calculate.py
python results-data/qwen-v2/merge_protocol.py
python results-data/three-dataset-tables/build_tables.py
python results-data/build_druglike_table.py
python results-data/ref-panel/build_supporting.py
python results-data/core-ref-protocol/build_energy_table.py
python plinder22-data/plot.py
python results-data/three-dataset-tables/plot_casf.py
python results-data/core-ref-protocol/build_figures.py
python results-data/size-flexibility/build_figure.py
python plinder22-data/molecule-tracking/plot.py
python build_candidate_figures.py
```

Copy the three `plinder22-molecule-tracking` artwork files from its source
folder to `figures/figure-plinder-molecule-tracking.{pdf,png,svg}` after rendering.
Do not rerun historical exporters to restore the former primary checkpoint.

Validation runs `qwen-v2/verify.py`, `ref-panel/verify.py`,
`plinder22-data/verify.py`, and `plinder22-data/qwen-dedup/verify.py`. They check
checkpoint-specific aggregates, exact probabilities, raw-coordinate energy
summaries, unchanged comparator measurements, table values and rank emphasis,
reference identities, molecule tracking, citations and figure links. The figure
artwork was inspected. The manuscript was not compiled locally, in accordance
with the paper workflow; compilation is reviewed through Overleaf.
