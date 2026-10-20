# PLINDER — cleaned per-ligand crystal conformers

```
plinder/
├── raw/        # TODO: original PLINDER download/version (not yet added)
├── processor/  # TODO: cleaning script + logs that produced cleaned/ (not yet added)
└── cleaned/
    ├── split_by_molecule_dedup_0.05_anchor_triplet.zip   # as supplied, 2026-10-20 (1,070 records, untouched)
    ├── split_by_molecule_dedup_0.05_anchor_triplet/      # working copy, 22 SDFs, 3 records removed (see below)
    └── stats/      # plinder_stats.py + its JSON output
```

Zip sha256: `32c4f40c83c83537553b9c2e66528d366e87bb7627a2aa944aa0827a198ee938`.
Provenance (PLINDER release, selection of the 22 ligands, cleaning code) is not yet recorded here.

## Contents

One SDF per ligand, named by CCD code. Every record is a heavy-atom-only (no H)
crystal ligand pose, RDKit-parseable, already deduplicated at 0.05 Å symmetric
RMSD and aligned on a fixed per-ligand anchor atom triplet (`out_align_mode=anchor_triplet`).

| | |
|---|---|
| ligands | 22 |
| conformers (deduplicated representatives) | 1,067 (8–117 per ligand, median 40); 1,070 as supplied |
| raw references before dedup | 2,261 (2,258 represented after the removal) |
| distinct PDB entries / PLINDER systems | 552 / 1,045 |
| heavy atoms / rotatable bonds / MW | 11–37 / 1–11 / 162–506 Da |

SD properties on every record: `molecule_name`, `ccd`, `reference_index`, `pdb_id`,
`system_id`, `duplicate_system_id_entry`, `dedup_cluster_id_rmsd0.05_symmetric`,
`dedup_cluster_size_rmsd0.05_symmetric`, `is_cluster_representative_rmsd0.05_symmetric`
(always True), `cluster_size`, `cluster_members_ref_idx`, `conf_index` (original IDs kept;
gaps at CLM 22, 23 and STO 37 after the removal),
`out_align_mode`, `anchor_triplet_atoms_0_based`.

Checked (2026-10-20, `stats/plinder_stats.py`):

- one heavy-atom graph per file; one anchor triplet per file; `ccd` matches filename;
- as supplied, `cluster_members_ref_idx` partitioned each ligand's raw references `0..N-1` exactly
  (after the removal, CLM 75/76 and STO 53 are unrepresented),
  each representative is in its own cluster, and `len(members) == cluster_size`;
- minimum pairwise in-place (no realignment) symmetric RMSD between kept
  conformers is ≥ 0.05 Å in every file, consistent with the dedup threshold.

Things to know:

- **3 wrong-stereo records removed on 2026-10-20** (likely PDB modelling errors). Stereo
  was perceived from the 3D coordinates; the affected centres were clearly tetrahedral
  (|signed volume| 2.0–2.5 Å³), so these were real inversions, not distortions.
  - CLM conf 22 and 23 (PDB `2uxp`, `2uxp__1__1.A__1.C` / `2uxp__1__1.B__1.D`, reference_index
    75 / 76): the (S,S) enantiomer; every other record is (R,R)-chloramphenicol.
  - STO conf 37 (PDB `2oic`, `2oic__1__1.A__1.E`, reference_index 53): one of four centres
    inverted (atom 3, R→S), an epimer.
  Each was its own cluster (`cluster_size` 1), so no other record's membership changed; those
  three reference indices are simply no longer represented. Removal was text-level: every
  remaining record is byte-identical to the zip. Each file now has a single stereoisomer.
  Per-record validity checks and symmetry-aware RMSD (graph matching ignores chirality) do
  not detect this; it needs a cross-record stereo comparison.
- **`system_id` is not a unique key:** 42 records (`duplicate_system_id_entry=True`) fall in
  20 systems containing 2–3 copies of the ligand; 1,070 records vs 1,048 systems. The copies
  differ by 0.13–5.2 Å RMSD (most > 1 Å, up to ~11 Å single-atom displacement), i.e. distinct
  ligand instances/poses in one system, as the cleaning notes state. Key on `(ccd, conf_index)`
  or `(ccd, reference_index)`. (Counts are as supplied; after the removal 1,067 records span
  1,045 systems.)
- **W11 (pleconaril) is dominated by symmetry copies:** 716 raw references (viral
  capsid copies) collapse to 8 representatives, one per PDB entry. The cleaning notes
  say 9 real experiments; the 9th is not present here (unresolved).
- Cleaning notes report 2,450 → 1,179 before removing bortezomib (covalent); this set's
  2,261 → 1,070 implies bortezomib contributed 189 raw / 109 kept records.

## Per-ligand statistics

`conf` kept conformers, `raw` references before dedup, `pdb`/`sys` distinct entries,
`hvy` heavy atoms, `rotb` rotatable bonds, `chiral` stereocentres. RMSD (Å, heavy-atom,
symmetry-aware) over all kept pairs: `minIP`/`medIP` in place (anchor-aligned frame),
`medB`/`maxB` after optimal superposition.

| ccd | name | conf | raw | pdb | sys | hvy | rotb | chiral | minIP | medIP | medB | maxB |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| STO | Staurosporine | 117 | 176 | 82 | 116 | 35 | 2 | 4 | 0.05 | 0.37 | 0.33 | 1.27 |
| BB2 | Actinonin | 82 | 99 | 38 | 79 | 27 | 11 | 3 | 0.13 | 1.20 | 0.84 | 3.11 |
| AZM | Acetazolamide | 81 | 92 | 37 | 81 | 13 | 2 | 0 | 0.06 | 0.60 | 0.38 | 1.57 |
| TCL | Triclosan | 81 | 151 | 32 | 79 | 17 | 2 | 0 | 0.07 | 0.35 | 0.21 | 1.93 |
| CLM | Chloramphenicol | 67 | 155 | 25 | 66 | 20 | 6 | 2 | 0.08 | 0.97 | 0.84 | 2.66 |
| CYZ | Cyclothiazide | 63 | 146 | 28 | 63 | 24 | 2 | 4 | 0.07 | 0.64 | 0.53 | 1.60 |
| NCT | Nicotine | 62 | 76 | 20 | 62 | 12 | 1 | 1 | 0.08 | 1.07 | 0.44 | 1.53 |
| TOP | Trimethoprim | 58 | 60 | 29 | 56 | 21 | 5 | 0 | 0.08 | 1.09 | 0.67 | 1.74 |
| DIF | Diclofenac | 50 | 55 | 20 | 45 | 19 | 4 | 0 | 0.10 | 1.47 | 1.02 | 1.64 |
| STI | Imatinib | 50 | 57 | 28 | 50 | 37 | 7 | 0 | 0.16 | 0.97 | 0.36 | 3.35 |
| EF2 | S-thalidomide | 40 | 73 | 23 | 40 | 19 | 1 | 1 | 0.08 | 0.43 | 0.18 | 0.47 |
| OHT | 4-Hydroxytamoxifen | 40 | 81 | 14 | 39 | 29 | 8 | 0 | 0.08 | 0.89 | 0.64 | 1.39 |
| IBM | IBMX | 38 | 57 | 20 | 38 | 16 | 2 | 0 | 0.10 | 0.99 | 0.72 | 1.17 |
| IMN | Indomethacin | 33 | 54 | 23 | 31 | 25 | 4 | 0 | 0.27 | 2.18 | 1.38 | 2.50 |
| TSN | Trichostatin A | 33 | 35 | 15 | 29 | 22 | 6 | 1 | 0.11 | 1.13 | 0.64 | 1.66 |
| DEX | Dexamethasone | 31 | 32 | 19 | 31 | 28 | 2 | 8 | 0.08 | 0.42 | 0.26 | 1.30 |
| 478 | Amprenavir | 29 | 30 | 25 | 29 | 35 | 11 | 3 | 0.14 | 0.50 | 0.31 | 2.75 |
| LOC | Colchicine | 29 | 30 | 19 | 29 | 29 | 5 | 1 | 0.10 | 0.41 | 0.30 | 0.79 |
| TPF | Fluconazole | 27 | 32 | 17 | 27 | 22 | 5 | 0 | 0.06 | 1.54 | 0.99 | 1.96 |
| FOM | Fosmidomycin | 26 | 28 | 15 | 26 | 11 | 5 | 0 | 0.08 | 1.06 | 0.65 | 1.38 |
| BMF | Bromosporine | 22 | 23 | 15 | 21 | 28 | 5 | 0 | 0.19 | 0.98 | 0.91 | 2.75 |
| W11 | Pleconaril | 8 | 716 | 8 | 8 | 27 | 6 | 0 | 0.56 | 1.80 | 0.93 | 1.24 |
