# Finalized PLINDER-22 minimized baselines

The manuscript adds `rdkit_random_minimized_22_druglike` and
`torsion_minimized_22_druglike` from the finalized 23 October 2026 runs. The
superseded fresh-draw, 23-molecule/UFF, and intermediate minimized runs are not
used. The display names are RDKit minimized and Torsion minimized.

## Sources and reproduction

The finalized panel CSVs and RMSD matrices are in
`/mnt/weka/vtarasov/code/casf-benchmark/data/results/`, under
`plinder_panel/<method>/` and `cache/druglike_rmsd_matrices/<method>.npz`.
Generation outputs and manifests are in
`/mnt/weka/vtarasov/outputs/casf_benchmark_runs/druglike_baselines/gen_results/`.
The existing eight cleaned-panel methods come from the mbedrosian benchmark
checkout, rather than the other checkout's older reference-panel exports.
All imported files and inspected calculation code have SHA-256 records in
`../provenance.json`. Generation and validity evaluator snapshots are retained
here along with both generation manifests. They are evidence snapshots, not
paper rendering entry points. Large generation pickles and full RMSD matrices
remain at their recorded source paths; `../nearest-distances.npz` archives all
nearest distances needed to reproduce the figures and coverage summaries.

Refresh the paper archive from the benchmark environment with:

```bash
python jcim-overleaf/plinder22-data/export.py \
  /mnt/weka/mbedrosian/code/casf-benchmark \
  --minimized-benchmark /mnt/weka/vtarasov/code/casf-benchmark \
  --baseline-runs /mnt/weka/vtarasov/outputs/casf_benchmark_runs/druglike_baselines
```

Then run the table and plot commands in the parent README, including
`python jcim-overleaf/results-data/size-flexibility/build_figure.py` for Figure 5.
`python jcim-overleaf/plinder22-data/verify.py` checks the ten-method archive,
tables, reference inventory, citations, includes, and unchanged CASF table data.
No new generation, optimization, alignment, or clustering is performed by the
paper export/rendering commands.

## What was optimized

Each new pool is a conformer-for-conformer minimization of its existing raw
counterpart, preserving order. Hydrogens are rebuilt with assigned coordinates,
MMFF94s runs for at most 500 iterations, and hydrogens are removed again. Both
converged and iteration-limit outputs are retained (status 0 or 1). No new
sampling, extra clash filter, or UFF fallback is used. Torsion retains its
experimental-reference seed. This is distinct from the separately sampled
raw/minimized CASF pools.

Both minimized output dictionaries contain 22 molecules and 22,000 conformers,
with 1,000 per molecule. The manifests' `n_molecules: 23` counts attempted
identities, including bortezomib with missing MMFF parameters and no output.
That molecule is already absent from the cleaned PLINDER reference panel.
Among the included molecules, 83 RDKit and 355 torsion conformers reached the
iteration limit; none failed or was dropped. These outputs remain in coverage,
precision, matching, and clustering, which use full supplied pools.

The RMSD cache shapes match the same 22 molecular identities and 1,067 cleaned
references. Per-molecule and aggregate COV/MAT values reproduce the finalized
exports within 1e-12. Full-pool cluster counts and occupancy reproduce their
summaries. Candidate identities and source hashes are recorded without
substituting a newer Qwen checkpoint.

## Corrected raw-baseline validity

The accompanying validity evaluation fixes atom ordering before comparison to
the reference topology. The previous raw-baseline records could report spurious
stereochemical failures. The corrected evaluator calls `renumber_to_reference`
before PoseBusters. The paper imports corrected raw PB counts and the associated
PB-subset energy records for the same 22 identities; the relevant selected rows
are archived here as `*-corrected-validity-energy.csv`.

The other checkout's raw geometric exports describe an older 23-molecule
reference collection and are **not** used for COV/MAT or clustering. Every
pre-existing cleaned PLINDER geometric measurement and all six learned-method
records are preserved. Corrected validity is measured separately and does not
filter the supplied pools. The raw-baseline PB percentages change in Table 4;
there is no claim that minimization universally improves validity.

## Displays and interpretation

Tables 1/2 fill the existing minimized-method PLINDER cells. Table 4 and S4 add
two rows, and their rankings are recalculated from unrounded values. Figure 1 keeps all ten methods in the top scatter panel, using the existing
palette and square markers for minimized classical variants. The bottom curve
panels show only RDKit raw and the six learned methods; RDKit raw is dashed.
Figure 5 fills the existing minimized-method PLINDER rows, leaving CASF cells
unchanged. Both torsion variants carry the reference-seeded marker.

Minimization reduces cluster counts and improves COV-R at 0.75 Å for both raw
pools, but decreases COV-R at 0.5 Å. Minimized torsion leads COV-P and MAT-P at the
primary cutoff, with its experimental seed explicitly qualified. Qwen still
leads primary-cutoff reference coverage and MAT-R; LoQI leads precision among
the learned methods. The six-learned-method review candidates retain their
stated subset and original measurements.
