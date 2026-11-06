# Big-data deduplication and CASF/PLINDER leakage filtering: methods description

Written for the paper's data and methods sections. Everything below describes what the code does, with counts taken from the evidence JSONs and run outputs.

- Part I covers the big-data corpus release: code in [coordtoken/data_collection/](../coordtoken/data_collection/).
- Part II covers the CASF-285 / PLINDER-23 leakage filter: code in [scripts/dedup_casf_plinder/](../scripts/dedup_casf_plinder/).

---

# Part I — Big-data corpus: splits, molecular separation, exact deduplication

**Final release:** `merged_deduplicated_20261006`, with **189,792,028 train / 1,802,906 val / 5,553,421 test** conformers. GEOM-Revisited (23,404 rows) is kept as a separate test file, which brings the test total to 5,576,825.

The stages run in this order: prepared snapshot → Group A → Group B → independent check → Group C → merge → exact-duplicate audit → exact deduplication → post-dedup audit and finalize.

## I.0 Molecular identity: one parser for every stage

All split decisions use one pinned parser, [audit_parser.py](../coordtoken/data_collection/group_a_split/audit_parser.py). Each stage checks the parser's SHA-256 and the RDKit version (2025.09.6) before it runs.

- **`inspect(text, vocab)`** tokenizes the enriched string with a strict regex. The token types are `[atom]<x,y,z>`, ring closures, bonds, branches, digits and `.`, and any unrecognised fragment raises an error.
  - It checks that each coordinate triplet has exactly 3 values, is finite, uses a standard numeric format, and fits in float32.
  - It removes CoordToken's decorative brackets from organic-subset atoms (`[C]`→`C`; the set is B C N O P S F Cl Br I b c n o p s) and rebuilds a plain SMILES string.
  - It requires the parsed atom count to equal the number of coordinates.
  - It records review flags that do **not** exclude a row: coincident atoms, bonds shorter than 0.7 Å or longer than 2.5 Å, out-of-vocabulary tokens, and InChI warnings.
- **`graph(smiles)`** parses with `sanitize=True, removeHs=False` and clears atom-map numbers. It then calls `rdinchi.MolToInchi` and accepts only return codes 0 or 1 with an `InChI=1S/` (standard) prefix, before converting with `InchiToInchiKey`.
- **Identity policy** ([split_policy.json](../coordtoken/data_collection/prepared_snapshot_20261005/split_policy.json)):
  - The identity is the full 27-character Standard InChIKey of the sanitized SMILES decoded from the stored string.
  - Encoded stereo, isotopes and charges are kept.
  - There is no salt stripping, no tautomer rule beyond what Standard InChI does itself, no stereo-free matching, and no scaffold matching.
  - All conformers of a molecule are kept.
  - A row whose key cannot be generated is **quarantined**. It is never treated as a new, non-overlapping molecule.

**Prepared input snapshot** ([prepared_snapshot_20261005/](../coordtoken/data_collection/prepared_snapshot_20261005/README.md)):

- The inputs are grouped CSVs. Some of them had already been filtered historically.
- ChEMBL3D was replaced by source-checked populations, so its membership differs from the defective backups.
- Ten PubChem rows (CID 141568402) were repaired from source and appended byte-for-byte.
- Nine rows are still quarantined, each identified by (source path, name, row SHA-256):
  - 6 ChEMBL rows (unresolved stereo normalization)
  - 1 PubChem row (graph parse failure)
  - 2 ZINC rows (1 parse failure, 1 with no InChIKey)
- The key of the quarantined ChEMBL **validation** row is still reserved as a validation identity.

**Ownership roles**, used in every stage: test = 3, val = 2, train = 1, unowned = 0. Registries merge with `max(role)`, so test takes precedence over val, and val over train.

## I.1 Group A: datasets with predefined splits

There are 14 datasets × 3 splits = 42 CSVs: BindingMOAD, BindingNet-High/Mid/Low, CrossDocked2020, DAVIS-3D, HiQBind, Kinodata-3D, OMol25-bio, PLINDER, SAIR, SPINDR, ChEMBL3D and GEOM.

- [prepare.py](../coordtoken/data_collection/group_a_split/prepare.py) pins the files, sizes, mtimes and hashes, plus the GEOM-Revisited, ∇²DFT-structures and ∇²DFT-scaffolds tests as fixed references. It splits each file into 64 MiB byte chunks.
- [process_group_a.py](../coordtoken/data_collection/group_a_split/process_group_a.py):
  1. **`verify_input`**: full SHA-256 of every input.
  2. **`extract`** (parallel per chunk): writes a row-aligned `keys.txt` sidecar, per-key counts, and failures. A key failure in a reference file (val/test) stops the run; failures in train are quarantined.
  3. **Registry build:**
     - Test key set = all Group A tests ∪ ∇²DFT structures ∪ ∇²DFT scaffolds ∪ GEOM-Revisited.
     - Val key set = all Group A val ∪ the reserved quarantined key.
     - The ∇²DFT **conformation** test is deliberately left out of the test set. Shared test-test keys are reported, not resolved.
  4. **`route`**:
     - Test rows are always kept.
     - A train or val row is excluded if its key is in the test set.
     - A train row is also excluded if its key is in the val set.
  5. **`write_chunk`**: streams source bytes against the key sidecar, writes kept and excluded rows with reasons, and checks that no retained key overlaps a higher-priority set. It also checks that no quarantined ChEMBL record reappears.
  6. **`assemble`**: concatenates the parts, re-reads the output to check its hash, and reconciles row counts.
- **Result:**
  - Train 11,094,343 → 10,339,711 (−754,632)
  - Val 901,495 → 874,019 (−27,476)
  - Test 765,424, unchanged

## I.2 Group B: datasets with no predefined split

KIBA-3D, OMol25-small, PubChem3D and ZINC are split **iteratively, in that order**.

- [prepare_index.py](../coordtoken/data_collection/group_b_split/prepare_index.py) and [index_inputs.py](../coordtoken/data_collection/group_b_split/index_inputs.py) key all four pools plus the ∇²DFT training pool: 183,260,646 rows with zero new failures. Chunk receipts are pinned by hash.
- [split_pools.py](../coordtoken/data_collection/group_b_split/split_pools.py):
  - **`seed_registry`** (SQLite `ownership.sqlite`) assigns roles:
    - role 3: Group A test keys
    - role 2: Group A val keys
    - role 1: Group A train keys and **∇²DFT train-pool keys**, so these can never become new Group B holdouts
  - **`allocate`** runs per dataset. It aggregates `key → n_rows`, looks up each key's owner, and sets the target T = round(0.005 × input rows).
  - **`choose`**, run for test first and then val:
    - Every key starts as train. Previously owned test keys start as `excluded_test_quota` and owned val keys as `excluded_val_quota`.
    - If the owned candidates (test, or val) fit within T, all of them are assigned.
    - Otherwise **`select_ranked`** orders whole key groups by `SHA-256(seed=42 ‖ dataset ‖ split ‖ key)` and takes groups until the row count reaches T. Overshoot is allowed because groups are never split.
    - Candidates that are not selected stay **excluded**; they never go back to training.
    - Any shortfall is filled only from **unowned** keys (role 0), using the same ranking. A shortfall that cannot be filled is reported, not forced.
  - After each dataset, its test/val/train keys are added to the registry, so later datasets inherit them.
  - **`write_chunk` / `assemble`**: byte-preserving routing, with row/key alignment, hash, and accounting checks.
  - **Final check:** every dataset's assignments are checked again against the *final* registry, which catches holdouts introduced by later datasets.
- **Result:**

| Dataset | Train | Val | Test | Excluded |
|---|---:|---:|---:|---:|
| KIBA-3D | 260,129 | 1,533 | 1,576 | 17,770 |
| OMol25-small | 4,033,245 | 21,597 | 37,384 | 226,716 |
| PubChem3D | 93,802,743 | 473,756 | 473,758 | 0 |
| ZINC | 74,401,034 | 375,763 | 375,763 | 0 |

- The OMol25-small test set overshoots (37,384 rows against a target of 21,595) because one InChIKey group alone has 24,255 conformers.

## I.3 Independent holdout check

[verify_holdout_precedence.py](../coordtoken/data_collection/group_b_split/verify_holdout_precedence.py) is a separate code path. It strips coordinates with a regex and computes keys with RDKit's direct `inchi.MolToInchiKey` (**`direct_key`**). It does not use the ownership registry.

- It re-parses all 5,707,138 holdout rows in 39 files.
- The test and val key unions share **zero** keys.
- Retained A/B training has zero matches to the union (checked through the saved key sidecars, not a fresh parse of all 183M rows).
- All 57 file hashes match.
- PubChem traces show all 270 keys that were in both the original test and val lists went to test.

## I.4 Group C: ∇²DFT

[prepare.py](../coordtoken/data_collection/group_c_split/prepare.py) pins the A/B evidence. [run_nabla.py](../coordtoken/data_collection/group_c_split/run_nabla.py) then:

1. Looks up each ∇²DFT key in the final Group B registry.
   - Role 3 → `excluded_test_overlap`.
   - Role 2 → existing-validation candidate.
2. **`blocked_by_b`**, plus a Group A scan, sets `other_train = 1` for keys found in any **retained** A or B training set.
3. **`select_validation`**:
   - Target = round(0.005 × (train pool + all three ∇²DFT tests)) = 63,325.
   - Existing-val matches are taken first, with ranked subsampling if they exceed the target.
   - The remainder is filled only from keys with `other_train = 0`. No other dataset's training rows are ever removed.
4. **`copy_test`**: the three ∇²DFT test files are copied byte-for-byte.
5. **Verification:**
   - The written val CSV is re-parsed and its keys must equal the selected set.
   - ∇²DFT train vs the actual A/B holdouts from I.3: 0 keys.
   - ∇²DFT val vs actual tests: 0 keys.
   - New ∇²DFT val vs A/B train: 0 keys.
   - All 54 A/B outputs are re-hashed.
- **Result:**
  - Train 8,689,383; 5,156 rows excluded for test overlap.
  - Val 63,340 = 29,684 existing-val matches + 33,656 new rows (15 over target).
  - The conformation test overlaps training **by design**: it evaluates new conformations of molecules seen in training.

## I.5 Merge with provenance

[merge.py](../coordtoken/data_collection/merge_release/merge.py) adds no new split decisions, shuffling, or deduplication.

- Output schema: `name, enriched_text, dataset, source_name, source_file_id, source_row_index, inchikey, test_subset`.
  - `name` = `<source_file_id>:<row_index>`.
  - `source_file_id` = `s` + the first 16 hex characters of SHA-256(path ‖ file sha).
- **`prepare`** plans one task per source chunk.
- **`run_chunk`** walks each source row and its kept/excluded split parts in parallel. Every source row must match exactly one part's next line, byte-for-byte. Hashes and counts are verified at every step.
  - Keys for the ∇²DFT conformation test, which was never indexed, are computed with the same frozen parser.
- **`copy_region`** pre-sizes each output file and `pwrite`s parts into disjoint offsets. **`verify_output`** re-reads each region and checks its hash and line count.
- **Result:** 191,526,245 / 1,810,008 / 5,560,978 rows.

## I.6 Exact-duplicate audit (read-only)

[audit.py](../coordtoken/data_collection/exact_duplicate_audit/audit.py):

- Equality means an exact match of the `enriched_text` field: no canonicalization, rounding, or name comparison. Different conformations stay distinct.
- **`scan_worker`** (64 workers): for each row it writes a record (SHA-256 of the text, category = split × dataset × subset, byte offset) into one of 256 buckets, chosen by the digest's first byte.
- **`analyze_bucket`** (32 reducers): sorts each bucket by digest, then fetches every candidate row from the files and confirms full-string equality. Hash collisions are counted separately and cannot produce a false duplicate.
- **Result:**
  - 198,920,635 rows scanned (including GEOM-Revisited).
  - 1,689,701 duplicate groups; 1,750,913 rows redundant globally.
  - Redundant within split: train 1,734,217, val 7,102, test 7,557.
  - **Zero** exact strings shared across splits.
  - GEOM-Revisited has no internal repeats, but 2,037 of its strings also appear in the GEOM test set.

## I.7 Exact deduplication, then finalize

[deduplicate.py](../coordtoken/data_collection/exact_deduplication/deduplicate.py):

- **`prepare`**
  - Refuses to run if the audit found any cross-split shared text.
  - For each duplicate group, within each merged file, it keeps the **lowest byte offset** (the first occurrence in file order) and records the others as drops, each with (offset, keeper offset, digest).
  - It checks that no offset is repeated, that no keeper is itself dropped, and that the drop counts equal the audit's counts.
- **`filter_region`**
  - For every dropped row it re-reads the keeper and requires the same full text, the same InChIKey, a matching digest, and a different record ID.
  - Each removed row is written to `duplicate_aliases.csv` as (removed_name, kept_name, original metadata, split).
- **`assemble`**: pwrite assembly, region readback, and whole-file SHA-256.

[finalize.py](../coordtoken/data_collection/exact_deduplication/finalize.py) runs only after a **fresh** audit of the deduplicated release. That audit must show:

- zero cross-split overlap,
- zero repeats within any dataset or within train/val,
- and that the only remaining matches are GEOM-Revisited ↔ GEOM test.

It then writes `release_status.json` and `SHA256SUMS`, and swaps the `merged_final` symlink atomically. [summarize_aliases.py](../coordtoken/data_collection/exact_deduplication/summarize_aliases.py) classifies each removal relative to its keeper:

- Train: 1,286,145 from the same source file; 448,072 from a different dataset (447,940 of those are OMol25-small → OMol25-bio).
- Val: 7,098 same file; 4 cross-dataset.
- Test: 7,554 same file; 3 cross-dataset.

**Per-dataset final counts** are in [paper_dataset_split_revision.md](../coordtoken/data_collection/paper_dataset_split_revision.md). That file already contains a methods draft for this part, and its numbers match the evidence.

**Scope to state in the paper:**

- Separation is molecule-level, by full InChIKey. It is not scaffold-level: scaffold-test protection excludes only exactly matching molecules.
- Deduplication removes exact strings only. It does not remove geometric near-duplicates.
- Group B quotas count input rows.
- The training CSV is in source order, not shuffled.
- The historical 190,214,248-row collection and the intermediate 189,925,494-row ∇²DFT revision in [FILTERING_METHODS.md](../coordtoken/data_collection/FILTERING_METHODS.md) are **superseded**. Don't cite them as final.

---

# Part II — CASF-285 / PLINDER-23 leakage filtering

The **input** is the deduplicated release above (train and val). The **queries** are:

- **CASF-2016 core**: 285 ligands from mol2 files (sanitized, H removed), listed in `ligand_descriptors.csv`.
- **PLINDER-23**: 23 ligands from the curated drug-like shortlist (`plinder_frequency/2024-06-v2/druglike_curated_shortlist.tsv`), parsed from SMILES.

There were four iterations. Each later one removes a superset of the previous one, through a "safety net". **The final training corpus is the K1/K2 ("protomer") variant**: `protomer/final/merged_train_clean_shuffled.csv`. It feeds `scripts/preprocess_casf_dedup` and the data key `casf_dedup_uniform_binned_pm12_train`.

## II.0 Shared machinery ([dedup_casf_plinder.py](../scripts/dedup_casf_plinder/dedup_casf_plinder.py))

**Identity**

- Every training row is decoded with **`decode_cartesian_v2`** (it rebuilds the SMILES and parses with sanitize on and Hs kept).
- **`canonical_smiles`** then removes Hs, calls `AssignStereochemistry(cleanIt, force)`, and writes canonical isomeric heavy-atom SMILES.
- **`identity`** re-parses that SMILES and computes the InChIKey and fingerprint from the re-parsed molecule. Queries and training rows go through exactly the same steps.

**Similarity** follows DODock:

- **`morgan_generator`**: 2048-bit, radius-2 Morgan fingerprint (ECFP4), chirality off, no count simulation.
- Tanimoto must be **strictly greater than 0.4**.
- It is computed on the molecule as written, with no charge or tautomer standardization.

**Supporting steps**

- **`prepare-queries`**: builds `queries.csv` with 308 rows. It checks that the SMILES round-trip leaves the fingerprint unchanged and notes any mismatch with the shortlist's InChIKeys.
- **Chunking**
  - **`chunk_bounds` / `iter_chunk_lines`** split the file into byte ranges (128 for train, 4 for val). A line belongs to the chunk that contains its first byte, so the scan and filter stages see identical rows in identical order.
  - **`parse_row`** has a fast path, falls back to the real CSV parser, and aborts the run on any unparseable row.
- **`scan`** → **`_scan_batch`**
  - Order-preserving multiprocessing (`imap`) with a per-batch SMILES cache.
  - Writes a row-aligned array of 27-byte keys (`keys.bin`).
  - For each distinct SMILES, computes Tanimoto to every query (`BulkTanimotoSimilarity`) and keeps the per-dataset maximum and its argmax when it exceeds 0.4.
  - Counts hits per source and logs mismatches against the CSV's `inchikey` column.
- **`merge`** (over train and val together):
  - Removal key set for each dataset = query keys ∪ keys of all hits.
  - These are reduced to **skeletons**: the first 14 characters of the InChIKey (the connectivity block). This excludes the stereo, isotopic and protonation layers.
  - CASF skeletons take precedence; PLINDER-only skeletons = PLINDER skeletons − CASF skeletons.
- **`filter_chunk`** sorts every row into one category:
  - 0 = kept
  - 1 = CASF-similar, moved to `casf_similar_<split>.csv` and **kept for analysis, not used for training**
  - 2 = PLINDER-only, **dropped**

  Rows are copied byte-for-byte.
- **`finalize`**: writes the CSVs, per-source tables, and the number of rows caught only by InChIKey (i.e., removed without being similarity hits themselves).
- **Shuffle**:
  - **`shuffle_scatter`** assigns each row to one of 256 buckets at random, using `default_rng([seed, task])`.
  - **`shuffle_plan`** computes byte offsets for each bucket.
  - **`shuffle_gather`** permutes each bucket in memory with `default_rng([seed, 10⁶+bucket])` and `pwrite`s it into place.
  - **`shuffle_verify`** checks bucket completion, total bytes, and newline count.
  - Random bucket assignment followed by a uniform in-bucket shuffle gives a uniform random permutation of all rows.
  - Seeds: 20261007 (plain), 20261009 (tautomer), 20261012 (final).

## II.1 Version history (useful for an ablation or footnote)

| Version | Match key | Train kept |
|---|---|---:|
| v1 | full 27-char InChIKey | 188,034,008 |
| plain | InChIKey[:14] | 187,927,881 (CASF 1,703,287; PLINDER-only 160,860) |
| tautomer ([tautomer_filter.py](../scripts/dedup_casf_plinder/tautomer_filter.py)) | InChIKey[:14] of `TautomerEnumerator.Canonicalize()`, in union with plain | — (used as a safety net) |
| **final K1/K2** | see II.2 | **187,398,101** |

- The tautomer variant's functions:
  - `tscan`: row-aligned tautomer keys plus a source-code byte per row.
  - `tmerge`: a funnel of exact / similarity / exact-of-similar matches, compared with the plain rule.
  - `casf-lookup`: matches against CASF core and REF ligands.
  - `tfilter --merge-dir`: applies the union rule.

## II.2 Final rule: K1/K2 normalized keys ([protomer_filter.py](../scripts/dedup_casf_plinder/protomer_filter.py))

**`normalized_keys`** works on a temporary copy of the molecule. The rows themselves are never modified.

- Normalization: `Uncharger(force=True, protonationOnly=True)` → `TautomerEnumerator.Canonicalize()`. If a step fails, it falls back to the molecule before that step, and the failure is counted.
- **K1** = InChIKey[:14] of the normalized molecule.
- **K2** = `rdMolHash` HetAtomProtomer of the normalized molecule with stereo removed, stored as a 16-byte BLAKE2b digest.
- Two molecules **match if K1 or K2 is equal.**

**Pipeline**

1. **`query-keys`**: K1/K2 for the 308 ligands.
2. **`pscan`** (train and val; same 128/4 chunking so arrays line up row-for-row with `keys.bin`): row-aligned K1, K2, and source code, plus every scanned molecule that exactly matches a ligand.
3. **`exact-set`**: builds **E** = the ligands plus their exact K1/K2 matches in the data. That is 285 + 675 for CASF and 23 + 100 for PLINDER. E is written out as a new query list.
4. **Similarity scan around E**: the `scan`/`merge` steps from II.0 run again, with E as the query list. This produced 381,812 similar molecules, up from 291,215 with the 308 ligands alone.
5. **`pmerge`** builds the removal sets for each dataset:
   - K1/K2 of the ligands, of E, and of every molecule with Tanimoto > 0.4 (original molecule) to any member of E.
   - The original full InChIKeys of those hits. These are matched directly as well, because normalization does not always reproduce a hit's K1.
   - Safety net: the plain and tautomer skeleton sets, so no row removed in an earlier version is returned to training.
   - Set sizes: CASF K1 235,351 / K2 235,606; PLINDER K1 28,437 / K2 28,435.
   - **`_pchunk`** computes the funnel and the transitions relative to earlier rules.
6. **`ref-exact`**: a fourth category for **CASF-2016 REF ligands with an exact ChEMBL3D intersection** (1,236 files; the `casf16_ref_chembl` evaluation set).
   - Only exact K1/K2 matches are removed; there is no similarity step for this set.
   - Matches go to their own file, `casf_ref_exact_train.csv`.
7. **`row_categories` → `pfilter`**: priority is CASF > PLINDER-only, then REF-exact applies only to rows that would otherwise be kept. Output is byte-preserving.
8. **`pfinalize`**: writes the CSVs and summaries. The rows are then shuffled with seed 20261012.
9. **`plookup`**: an audit of which CASF core/REF complexes (5,123 mol2 files) have a removed molecule: 1,754 complexes, 275 of them in the core set.

**Train funnel (CASF)**

| Step | Molecules | Rows | Notes |
|---|---:|---:|---|
| 1. Exact match to a ligand | 734 | 197,631 | covers 275 of 285 ligands |
| 2. Similarity to E | 315,706 | 1,944,085 | 80,299 molecules / 539,530 rows found only because of the expansion to E |
| 3. Exact match of a similar molecule | 9,693 | 59,940 | |
| Safety net only | | 7 | |

**Train funnel (PLINDER)**

| Step | Molecules | Rows | Notes |
|---|---:|---:|---|
| 1. Exact match to a ligand | 105 | 2,169 | covers 19 of 23 ligands |
| 2. Similarity to E | 42,916 | 167,699 | |
| 3. Exact match of a similar molecule | 2,090 | 4,075 | |

**Final train:** 189,792,028 →

| Category | Rows |
|---|---:|
| Kept | **187,398,101** |
| CASF-similar | 2,201,663 |
| PLINDER-only (dropped) | 163,118 |
| CASF-REF exact | 29,146 (1,484 molecules; 858 of 1,236 REF complexes) |

- Total removed: 1.26%.
- Compared with the plain rule, 497,556 rows moved from kept to CASF and 3,078 from kept to PLINDER-only. **No row moved from removed back to kept.**

**Per-source removal (final train):**

| Source | Rows | Kept | CASF-similar | PLINDER-only | CASF-REF exact | Removed % |
|---|---:|---:|---:|---:|---:|---:|
| BindingMoad | 20,480 | 15,157 | 4,192 | 180 | 951 | 25.99 |
| BindingNet-High | 210,402 | 197,256 | 11,192 | 1,293 | 661 | 6.25 |
| BindingNet-Low | 251,400 | 241,662 | 8,282 | 1,038 | 418 | 3.87 |
| BindingNet-Mid | 144,014 | 137,778 | 5,210 | 710 | 316 | 4.33 |
| CrossDocked2020 | 78,233 | 58,472 | 15,776 | 596 | 3,389 | 25.26 |
| DAVIS-3D | 129 | 56 | 55 | 18 | 0 | 56.59 |
| HiQBind | 24,495 | 18,347 | 4,412 | 314 | 1,422 | 25.10 |
| KIBA-3D | 260,129 | 247,088 | 11,497 | 342 | 1,202 | 5.01 |
| Kinodata-3D | 61,112 | 58,484 | 2,360 | 161 | 107 | 4.30 |
| OMol25_bio_mols | 1,314,710 | 1,106,575 | 198,091 | 1,797 | 8,247 | 15.83 |
| OMol25_small_mols | 3,551,242 | 2,777,850 | 767,897 | 1,200 | 4,295 | 21.78 |
| Plinder | 122,072 | 94,914 | 23,414 | 1,219 | 2,525 | 22.25 |
| SAIR | 1,521,900 | 1,485,057 | 30,049 | 5,963 | 831 | 2.42 |
| SPINDR | 25,340 | 20,374 | 3,938 | 320 | 708 | 19.60 |
| chembl3d | 1,428,111 | 1,398,651 | 24,416 | 4,480 | 564 | 2.06 |
| geom | 5,097,261 | 5,016,456 | 72,623 | 6,020 | 2,162 | 1.59 |
| nablaDFT | 8,689,383 | 8,623,394 | 62,978 | 2,968 | 43 | 0.76 |
| pubchem3d | 93,724,490 | 92,818,460 | 785,778 | 118,994 | 1,258 | 0.97 |
| zinc | 73,267,125 | 73,082,070 | 169,503 | 15,505 | 47 | 0.25 |

## II.3 Points to state precisely in the paper

1. **Val/test coverage.**
   - The K1/K2 rule was applied to **train only**.
   - Val was cleaned with the plain InChIKey[:14] rule: 1,802,906 → 1,779,945 kept (20,910 CASF-similar, 2,051 PLINDER-only dropped). The CASF-dedup training configs don't use a validation set (`dataset_path = ""`).
   - Test splits were not CASF-filtered.
   - The removal sets were built from hits in train **and** val.
2. **Expansion depth is bounded.** The chain is ligand → exact matches (E) → similar to E → exact matches of those. Similarity is not followed transitively beyond that.
3. **Two identity pipelines differ slightly.** The CASF stage's key (RemoveHs + forced stereo assignment + re-parse) disagrees with the release's `inchikey` column on 1,812 train rows and 9 val rows.
4. **Decode failures** would leave a row with an empty key that never matches, so it would be kept. The logs record zero decode failures and zero missing keys, so this didn't happen.

---

## Draft methods paragraph (CASF/PLINDER)

> To prevent ligand leakage into the CASF-2016 core set (285 ligands) and our PLINDER-23 benchmark, we filtered the deduplicated training corpus with a three-step rule. Each molecule was decoded from its coordinate-enriched string and assigned two normalized structure keys, computed on a copy after neutralization (RDKit Uncharger, protonation only) and tautomer canonicalization: K1, the connectivity block of its standard InChIKey, and K2, the stereo-free RDKit HetAtomProtomer hash. Two molecules match if either key is equal. We first collected all training molecules matching any benchmark ligand (set E). We then flagged every molecule whose ECFP4 fingerprint (2048 bits, radius 2, computed on the unmodified molecule) had Tanimoto similarity > 0.4 to any member of E, following DODock, and finally all molecules matching a flagged molecule by K1 or K2. Every conformer of a flagged molecule was removed. We also removed exact K1/K2 matches to the 1,236 CASF-2016 refined-set ligands used in our CASF-REF evaluation. This removed 2,393,927 of 189,792,028 training conformers (1.26%), leaving 187,398,101, which were shuffled with a fixed seed.
