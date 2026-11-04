# Independent PLINDER torsion update

The raw and minimized PLINDER torsion results now use:

- `torsion_independent_raw_druglike`
- `torsion_independent_minimized_druglike`

These replace `torsion_raw_druglike` and `torsion_minimized_22_druglike`, whose
starting coordinates came from an experimental reference. The replacement
removes the dagger and reference-seed qualification from the current manuscript.
Historical exports and dated provenance remain as history.

## Executed procedure

Each molecular graph is built from isomeric SMILES with no coordinates. Five
ETKDGv3 starts produce 200 accepted torsion perturbations each, giving 1,000 raw
conformers per molecule. Torsions are displaced within ±120 degrees, with severe
clash rejection. Experimental coordinates are not used; a coordinate-free
reference graph is checked only for molecular identity.

The minimized run preserves raw-conformer order and optimizes each member with
MMFF94s after adding hydrogens, for up to 500 iterations. Stereo-inverted results
are retried with chirality-preserving constraints. Converged and iteration-limit
results are retained after identity checks. All 22 molecules have 1,000 raw and
1,000 minimized outputs. `seed-provenance.json` retains per-start seeds, attempts,
counts, minimization-status summaries and output hashes; source hashes identify
the full upstream per-conformer records. Generation code and summary manifests
are archived here, but generation and minimization were not rerun for this edit.

## Results and scope

| Metric | Raw, previous → current | Minimized, previous → current |
|---|---:|---:|
| COV-R at 0.75 Å (%) | 83.0 → 80.1 | 88.2 → 84.4 |
| COV-P at 0.75 Å (%) | 51.1 → 36.2 | 65.6 → 56.7 |
| MAT-R (Å) | 0.461 → 0.532 | 0.484 → 0.494 |
| MAT-P (Å) | 0.891 → 1.053 | 0.717 → 0.876 |
| Mean clusters, 0.5 Å | 366.7 → 424.2 | 113.8 → 138.6 |
| Mean clusters, 1.0 Å | 135.7 → 153.6 | 52.3 → 69.1 |

Tables 1, 2, 4 and S4 and Figures 1 and 6 use the replacement results. Rankings
are recomputed from unrounded values; LoQI now leads COV-P and MAT-P, while
minimized torsion is second for COV-P and minimized RDKit is second for MAT-P.
Raw torsion slightly exceeds FlowR at the finer clustering radius; FlowR still
leads at 1.0 Å. All Qwen comparisons and other methods are preserved.

The source panel has the same 22 identities and 1,067 reference conformations.
Matching and geometric diversity use full supplied pools; PB remains separate.
Cached reference/generated RMSD minima reproduce all published COV/MAT values
and curves within 1e-12. No alignment or energy calculation was rerun.

The separate CASF seed audit reports 94 core and 1,236 original ref entries,
all with stored ChEMBL3D topology seeds and status `unchanged`. Its frozen
1,044-entry ref subset is unchanged. `casf-seed-audit.csv` normalizes escaped
line separators in the supplied manifest; its original hash is in the provenance.
The PLINDER replacement does not alter CASF tables, energy or sampling data.

## Reproduction

From `jcim-overleaf/`, with the benchmark environment activated:

```bash
python plinder22-data/independent-torsion/update.py /path/to/casf-benchmark
python results-data/three-dataset-tables/build_tables.py
python results-data/build_druglike_table.py
python plinder22-data/plot.py
python plinder22-data/molecule-tracking/plot.py
```

Copy the three `plinder22-molecule-tracking.{pdf,png,svg}` files from
`plinder22-data/molecule-tracking/` to
`figures/figure-plinder-molecule-tracking.{pdf,png,svg}` after rendering, and
refresh the rendered-file checksum manifests. `update.py` preserves every
unaffected CSV observation and method identity. It does not rerun scientific scoring.

Validation includes `verify.py` here, `plinder22-data/verify.py`, and
`plinder22-data/qwen-dedup/verify.py`. These check unchanged comparator and CASF
measurements, cached minima, table values/rank emphasis, figure records, citations,
cross-references and included artwork. The affected figures were visually inspected.
The full manuscript is reviewed in Overleaf and was not compiled locally.
