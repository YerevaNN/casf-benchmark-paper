# Common-hydrogen energy comparison

## Primary Qwen checkpoint update, 4 November 2026

All current primary Qwen results use `qwen_1p7b_v2_coordtoken_step41054`, displayed
as **Qwen 1.7B FSQ**. It replaces step47023 in core94, ref1044 and PLINDER22;
the separate Qwen 0.6B CASF-dedup result remains. Tables 1–4, S1–S4/S8/S9,
Figures 1–6/B1 and the optional review figures use the new values. Rankings
exclude the stored ChEMBL3D baseline as before. LoQI now leads PLINDER COV-R
at 0.75 Å; Qwen retains the lowest MAT-R. Training exclusions still require
an executed manifest and must not be inferred from the run name.

Energy and sampling were not published for this checkpoint, so they were measured
from its saved PB-passing conformers under the existing paper protocols. No
coordinate optimization or validity replay was performed. The current energy
totals are 78,058/78,058 core and 734,951/734,952 ref; the one undefined ref
energy is still FlowR's. Ref has 43 empty method–entry energy pools.

See `results-data/qwen-v2/README.md` for exact source identities, reproduction
and validation. This update supersedes older checkpoint statements below;
dated records retain their original lineage. Do not restore old Qwen values
by rerunning a historical exporter.


This analysis supplies the energy Results figure at the ChEMBL-count target on
core94. It uses the exact Qwen 1.7B FSQ step47023 ensemble and the same ten
pipelines plus ChEMBL3D-PB as Tables 1–2. It is a new energy calculation on the
archived retained ensembles, not a new generation run or a benchmark-excluded
Qwen evaluation.

The archived energy routine evaluates each molecule as stored, reusing a saved
post-minimization energy for marked minimized records. Inspection found that
Qwen output lacks explicit hydrogens while other generators generally include
them. Those archived values therefore do not supply a consistent comparison.
They are retained in the audit CSV for traceability, not used in the new figure.

`rescore.py` removes explicit hydrogens and reconstructs their coordinates from
each original heavy-atom structure using RDKit. It then evaluates a single-point
MMFF94s energy without optimizing any coordinates. Every method uses this same
hydrogen placement. Heavy-atom coordinates are checked for zero displacement;
previously stored energy properties are ignored. No molecule files or databases
are modified. The minimized baseline ensembles retain the geometries produced by
their generation pipelines, but receive no further minimization during scoring.
RDKit's [force-field documentation](https://www.rdkit.org/docs/source/rdkit.ForceField.rdForceField.html)
defines the energy calculation and kcal/mol units.

This protocol supersedes the earlier hydrogen-relaxed calculation at the author's
request. That calculation is preserved in paper commit `55ee5fa`; its values must
not be mixed with the current unoptimized results. The original archived energies
also remain unsuitable for this comparison because hydrogen representation varied.
The current calculation measures raw heavy-atom geometries with a common,
unoptimized hydrogen placement, not the original hydrogen coordinates supplied
by each generator.

The 94 stored ChEMBL3D reference ensembles are reconstructed using the repository's
identity-matched loader. All 94 have zero PoseBusters rejections in the pinned
snapshot. Counts and archived median/SD fingerprints are checked before rescoring.
Generated ensemble counts must match the database. Empty outputs are retained in
the audit; missing energy is not replaced by zero.

The within-ensemble standard deviation uses the population definition (ddof=0).
The plot summarizes these deviations by their median across molecules, alongside
mean 0.5 Å cluster counts for the same entries. The supporting table also reports
the mean of the energy SDs, median per-molecule mean energy, and median paired
difference of mean energy relative to the same molecule's ChEMBL3D ensemble.
This paired difference is not obtained by subtracting two aggregate medians.

A small SD measures narrow energy dispersion; it does not establish low absolute
energy or the absence of extreme conformers. Taking a median across molecules
limits the influence of exceptional molecules on the displayed summary, not on
the underlying dataset. These are force-field potential energies, not binding
free energies or equilibrium conformational populations.

Run in the analysis checkout:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python docs/energy_analysis/rescore_fast.py
python docs/energy_analysis/build_figure.py
```

The script records source-file digests, database and mapping digests, RDKit
version, protocol, calculation status, and maximum heavy-atom displacement.
`per_conformer.csv.gz`, `per_molecule.csv`, and `summary.csv` preserve all three
aggregation levels. The cache is an execution aid and is not needed to interpret
the archived results. The manuscript copies of these scripts retain the original
analysis-checkout path assumptions.

`rescore_fast.py` uses batched Zarr reads for the reference ensembles, preserving
the original candidate ordering, stereoisomer and atomic-number checks, and
reference sampling. It avoids repeatedly decompressing large coordinate chunks.
The original loader and optimized loader share the same energy/count checks;
validated cached ensembles are reused without changing their calculations.
