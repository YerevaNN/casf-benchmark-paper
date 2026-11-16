# PLINDER energy addition to Figure 2

Saved per-molecule energy SDs are aggregated by their median across molecules.
No energies, geometries, PoseBusters membership or clustering are recalculated.
The archived evaluator matches the recorded source hash; its MMFF94s scoring and
energy-summary functions match both current evaluation checkouts. `export.py`
checks the imported energy records against upstream records for all ten methods.
The protocol removes and rebuilds hydrogens, evaluates a single point without
optimization, uses nonbonded threshold 100 and ignores interfragment interactions.
Population standard deviation uses ddof=0.

Energies use PB-passing subsets; PLINDER cluster counts use full supplied pools.
Both summaries use molecules with defined energy SD: 21 for FlowR (no passing
4-hydroxytamoxifen conformer), 22 for every other method. This is why FlowR's mean
cluster count here differs from Table 1. All 192,356 passing conformers have finite
energy. Individual conformer energies and per-molecule means were not exported by
this PLINDER pipeline, so mean-energy and ChEMBL paired-difference columns in S8
remain CASF-only. The per-molecule SD summaries suffice to reproduce the figure.

```bash
python jcim-overleaf/results-data/plinder-energy/export.py
python jcim-overleaf/results-data/plinder-energy/render.py
python jcim-overleaf/results-data/core-ref-protocol/build_energy_table.py
python jcim-overleaf/results-data/plinder-energy/verify.py
```

Only export requires access to the recorded evaluation checkouts. Rendering uses
archived CSVs. Figure 2 has a logarithmic energy axis for PLINDER because raw
torsion spans a much wider range. CASF values remain unchanged. Figure 3 is unchanged.
Protocol/source checks are in `protocol-audit.json`; rendering input/output hashes
are in `rendered-files.json`. Source snapshots are evidence, not pipeline entry points.
