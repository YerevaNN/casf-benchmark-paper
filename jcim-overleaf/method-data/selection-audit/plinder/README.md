# PLINDER identity-selection evidence, 6 November 2026

The author description is preserved in `../../plinder23-selection-provenance.md`.
This archive adds the original frozen shortlist and independently checks its
frequencies and the current 22 identities against the locally available source tables.
`sources.json` records original and archived hashes; `checks.json` records the checks.
The large source tables are identified by hash rather than copied here.

The source table has 690,692 ligand-instance rows and 51,280 complete InChIKey groups.
The saved physicochemical screen has 25,392 rows. All three saved frequency fields
for all 23 shortlisted identities agree with direct aggregation of the source table.
All 22 evaluated identities and SMILES match the shortlist; bortezomib is the sole
removed identity. No physicochemical screen or experimental-coordinate cleaning was
rerun. Current tables and scores are unchanged.

## September filtering workflow is separate from v2 training provenance

The local September 14 copies of `dodock_similarity_filter_audit.py` and
`filter_serialized_jsonl_by_smiles.py` are archived here with `queries.csv`.
The latter contains exactly 285 CASF and 23 PLINDER queries, with the PLINDER queries
matching canonicalized shortlist SMILES. No PoseBench or Platinum queries are present.

The audit script finds exclusions with 2,048-bit radius-2 Morgan fingerprints,
chirality disabled, and Tanimoto strictly greater than 0.4. It supports FPSim2 and
NumKong search and records closest queries, similarity, source attribution and counts.
The separate writer loads `geom_smiles`, compares it with JSONL `canonical_smiles`,
removes every matching record, keeps retained records byte-for-byte, rebuilds uint64
byte-offset indices and writes counts. Its resume mode skips completed output files.
Finding excluded molecules and writing a filtered corpus are distinct steps.

These scripts do not establish the exact training input of a particular checkpoint.
The four evaluated v2 models instead have executed October corpus and normalized-key
exclusion manifests in `../../training-audit/`. Do not substitute this simpler
September procedure for that documented pipeline or use it to infer the older
step-64,000 checkpoint's lineage. That mapping remains a TODO.

The remaining PLINDER construction TODO concerns extraction of coordinates, their
release, protonation/tautomer/alternate-location handling, and the exact 0.05 Å
representative-selection procedure. Identity-selection provenance does not resolve it.
