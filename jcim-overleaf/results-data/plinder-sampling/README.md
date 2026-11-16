# Figure 5: PLINDER and CASF sampling

The main figure contains PLINDER COV-R at 0.75 Å and MAT-R above CASF core/ref
recovery at 0.75 Å. All four methods use the current paper identities, including
Qwen 1.7B FSQ v2 CoordToken step41054. No geometries were generated or optimized.

Render from archived curves, without the benchmark datasets:

```bash
python jcim-overleaf/results-data/plinder-sampling/render.py
python jcim-overleaf/results-data/plinder-sampling/verify.py
```

`calculate.py` optionally reproduces the PLINDER curves in the benchmark environment
with NumPy, pandas, SciPy, and `casf_benchmark` installed. It reads the original RMSD
cache paths in `provenance.json` and requires their hashes to match. Large matrices
remain in those source caches. No RMSD alignment is repeated. The archive contains
per-reference hit counts, per-molecule curves, and molecule-mean curves. MAT-R uses
the exact minimum-rank probabilities for a uniform subset without replacement;
COV-R uses the exact probability of at least one qualifying conformer. Expectations
are averaged within molecule before averaging equally across 22 molecules.

Budgets use the same 63-point grid (10–1,000) as CASF. Each budget is capped at
its available pool. PLINDER uses supplied conformers without a common PB gate;
CASF includes rejected candidates in its budget. These are sampling-efficiency
curves within evaluated pools, not generation-time or compute-efficiency curves.
The five nonfinite pairs have at most one per reference, so every displayed
subset necessarily has a finite match. No reference is excluded.

Validation checks monotonicity, 88 per-molecule full-pool endpoints, 4,268 reference
hit-count records, equal molecule weighting, current source identities and figure
hashes. The CASF curves are reused unchanged. `rendered-files.json` records all
rendering inputs and PDF/PNG/SVG output hashes. Methods and the evaluation appendix
document the calculation without formulas.
