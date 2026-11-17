# Protein-context pharmacophore-query analysis

This archive supports the pharmacophore-query Methods, Results, and Supporting
Figure. It contains compact summary tables and frozen provenance, not the large
candidate-level membership matrices or RMSD caches.

The PLINDER-22 pilot uses 932 eligible positional queries from 22 molecules and
1,067 cleaned reference complexes. The algorithmically selected PLINDER-ATC48
extension contains 48 molecules and 379 complexes; 311 queries from 44 molecules
are eligible. The extension broadens structural coverage but is not an unseen-
molecule benchmark: its exact isomeric InChIKey audit finds Qwen training records
for 17 of the 44 eligible identities and does not exclude near-neighbor or scaffold
overlap.

Render the archived summaries without loading conformers or recalculating matches:

```bash
python jcim-overleaf/results-data/pharmacophore-query/build_overlap_summary.py
python jcim-overleaf/results-data/pharmacophore-query/render.py
python jcim-overleaf/results-data/pharmacophore-query/verify.py
```

The primary endpoint uses 0.5 Å positional spheres and a requested budget of
1,000 accepted cached conformers. Exact uniform-subset probabilities are averaged
over queries within molecule and then equally across molecules. Confidence
intervals resample molecules. The 50:50 library allocates 500 draws to Qwen and
500 to LoQI. PoseBusters-gated summaries retain failed candidates in the sampling
denominator and do not replenish them.

The extension Qwen source is the accepted merged pool (47,943 conformers;
992--1,000 per molecule), so requested single-library budgets are capped at
availability. LoQI supplies exactly 1,000 conformers per molecule. Its CPU run
uses four deterministic 250-conformer shards for 38 molecules and eight
125-conformer shards for 10 slower molecules; the archived generation summary
and provenance record preserve the final layout and pool hash.

`provenance.json` records source paths, hashes, configuration identities, matching
completeness, and the training-audit scope. Candidate membership, native-complex
maps, conformer pools, PoseBusters masks, and verified RMSD matrices remain in the
benchmark repository paths identified there.

The `pilot-*` and `atc48-*` files are frozen copies of compact analysis outputs.
`build_overlap_summary.py` derives the exact-training-identity strata from the
archived per-molecule curves, and `render.py` writes the PDF, editable SVG, and PNG
artwork plus their input/output hashes. `build_provenance.py` is the collection-time
check used to confirm that each frozen copy is byte-identical to its recorded
benchmark source; it is not required to render or verify the self-contained archive.
