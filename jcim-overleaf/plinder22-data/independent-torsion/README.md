# PLINDER torsion results and procedure

Current source labels:

- `torsion_independent_raw_druglike`
- `torsion_independent_minimized_druglike`

## Executed procedure

Each molecular graph is built from isomeric SMILES. One ETKDGv3 starting
conformer supplies 1,000 accepted raw torsion perturbations. All eligible
torsions receive uniform displacements within ±120 degrees. The severe-clash
check uses 0.7 times the RDKit distance-geometry lower bound for heavy-atom
pairs, excluding bonded and angle-connected pairs. Experimental coordinates
are not used. The global seed is 42; per-molecule embedding and separate raw
and minimized sampling seeds are recorded in the source and manifests.

The minimized ensemble is a fresh torsion sample from the same ETKDG geometry,
not a transformation of the saved raw ensemble. It collects 1,500 candidates
passing the clash check, applies MMFF94s for up to 500 iterations, and repeats
the same clash check. Stereo mismatches are rejected. Additional tranches of
500 candidates fill the pool if necessary. Converged and iteration-limit
outputs are eligible; other minimization failures are rejected. The active
path does not use the chirality-constrained retry helper retained in the code.

All 22 molecules have 1,000 raw and 1,000 minimized outputs. The records in
`seed-provenance.json` retain seeds, attempts, rejection counts and output hashes.
The executed protocol is `fresh_torsion_mmff94s_clash_0.7`. The archived code
includes its shared generation/minimization implementation. Generation,
optimization and RMSD alignment were not rerun for the manuscript update.
PLINDER RDKit remains a paired raw-to-minimized comparison.

## Current results and display scope

| Metric | Raw | Minimized |
|---|---:|---:|
| COV-R at 0.75 Å (%) | 66.3 | 78.0 |
| COV-P at 0.75 Å (%) | 37.7 | 55.0 |
| MAT-R (Å) | 0.654 | 0.562 |
| MAT-P (Å) | 1.020 | 0.847 |
| Mean clusters, 0.5 Å | 359.0 | 122.3 |
| Mean clusters, 1.0 Å | 133.0 | 62.1 |

Tables 1, 2 and 4, main Figures 1 and 7, and Supporting Figures S5 and S6
use these measurements. Table S5 and the evaluation appendix describe the
procedure. Rankings use unrounded values, excluding the ChEMBL3D baseline.
FlowR leads PLINDER clustering at both radii; minimized RDKit is second for
COV-P and minimized torsion is second for MAT-P. The molecule-tracking figure
is ordered by mean MAT-R and retains its approved nine-method scope.

The panel contains 22 identities and 1,067 references. Matching and geometric
diversity use full supplied pools; PB remains separate. Cached RMSD minima
reproduce the published COV/MAT values and threshold curves within 1e-12.
All other methods' measurements remain unchanged.

`published-scope.json` records independent checks of the current staged public
release, `plinder-torsion-minimized-v1`, against its CASF source release. The
CASF dashboard is byte-identical; all 21 non-PLINDER/non-druglike extended
SQLite tables have identical rows. The CASF master CSV, per-ligand CSV and
base dashboard also match the recorded pre-generation hashes. Both public
PLINDER torsion rows agree with their source summary, per-molecule, recovery
and sampling-budget CSVs. Core94 and the paper's frozen ref1044 are unchanged.
The seed audit covers core94 and the 1,236-entry ref collection from which
that subset was selected, all using stored ChEMBL3D topology coordinates.

## Reproduction and validation

From the paper repository, with the benchmark environment activated:

```bash
python jcim-overleaf/plinder22-data/independent-torsion/check_published_scope.py /path/to/casf-benchmark
python jcim-overleaf/plinder22-data/independent-torsion/update.py /path/to/casf-benchmark
python jcim-overleaf/results-data/three-dataset-tables/build_tables.py
python jcim-overleaf/results-data/build_druglike_table.py
python jcim-overleaf/plinder22-data/plot.py
python jcim-overleaf/plinder22-data/molecule-tracking/plot.py
python jcim-overleaf/appendix-data/build.py
```

Copy `plinder22-molecule-tracking.{pdf,png,svg}` from the tracking directory
to `jcim-overleaf/figures/figure-plinder-molecule-tracking.{pdf,png,svg}`.
Refresh hashes for changed inputs/artwork after verifying that unaffected
measurements and displays remain identical. The landscape figure uses only
unchanged LoQI, FlowR and Qwen distances; its shared archive hash changes,
but its source entries and rendered artwork do not.

Run `independent-torsion/verify.py`, `plinder22-data/verify.py`,
`plinder22-data/qwen-dedup/verify.py`, `appendix-data/verify.py`, and
`results-data/conformer-landscapes/verify.py`. These check values, ranks,
unchanged comparators, included figures, citations and references. The full
manuscript is reviewed in Overleaf; no local manuscript compilation is needed.
