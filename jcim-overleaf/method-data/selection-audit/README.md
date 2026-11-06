# Verified panel-selection counts

The two unmodified matching logs record the executed core/ref selection and excluded
ligand identifiers. `sources.json` preserves their source locations and hashes.
The matching implementation is `src/casf_benchmark/cli/match_casf16_chembl3d.py`
in the benchmark repository, with stereochemical checks in `chembl3d/identity.py`.
The log category `missing_topology_sdf` also covers failure to find a topology
record with matching stereochemistry. `no_rotatable_bonds` tests eligibility
on the prepared topology; the paper uses wording that includes preparation failure.

Core: 285 - 0 - 189 - 1 - 1 = 94.
Ref: 4,838 - 11 - 3,566 - 12 - 13 = 1,236; then 1,236 - 192 = 1,044.
The final selection is recorded in `results-data/ref-panel/selection.csv` and
`selection.json`. No matching, generation, scoring, or filtering was rerun.

The source and preparation of the ref input archive remain unresolved.
Selected entry 6EL5 was deposited in 2017 and released in 2018, so the local
CASF16_REF name cannot identify a 2016 source release:
https://www.rcsb.org/structure/6EL5
No release attribution has been inferred from this observation.

PLINDER selection and cleaning evidence remains in
`method-data/plinder23-selection-provenance.md` and
`plinder22-data/reference-source-notes.md`. Upstream cleaning code and its
protonation, tautomer, and alternate-location rules remain explicit TODOs.
