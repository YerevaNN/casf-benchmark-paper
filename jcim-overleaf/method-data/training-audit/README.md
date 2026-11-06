# Training and evaluation evidence, 6 November 2026

The scientific account is split deliberately:

- **Main Methods (`../../methods.tex`)**: what CoordToken represents, the two training stages, evaluated Qwen identities, benchmark exclusion scope, basic sampling settings, and the qualification on representation comparisons.
- **Training appendix (`../../training-appendix.tex`, Table S8)**: architecture, checkpoint-specific training, exact exclusion rules/counts, corpus composition, Qwen optimizer and completed exposure, seeds, and outstanding records.
- **Evaluation appendix (`../../evaluation-protocol-appendix.tex`)**: validity configuration, alignment implementation/fallback, and molecule-level display calculations.

Both appendices are included in the main review PDF and standalone SI. No generation,
optimization, alignment, or scientific scoring was rerun. Main result numbers and
figure artwork are unchanged.

## Evidence hierarchy

The author-supplied *CoordToken: Tokenization of 3D Molecules for Autoregressive
Generation and Geometry Refinement* describes the architecture and the wider study.
Its full-corpus checkpoint and headline reconstruction results are **not** the
filtered checkpoint used here. The supplied PDF is identified by its hash in
`sources.json`; it is not redistributed in this archive. Its title, authors, year,
and workshop are transcribed from the supplied PDF. The final public bibliographic
record remains an author follow-up. Qwen3 report metadata were checked at
https://arxiv.org/abs/2505.09388 on 6 November 2026.

The filtered checkpoint is identified by its release, both epoch manifests,
completed results, input manifest, and epoch schedules. The first and second epoch
share the input-manifest hash. `coordtoken-input-manifest.json` points to the exact
retained CSV and its hash, rather than relying on matching record counts.
`fsq-encoding-manifest.json` binds the Qwen export to this same source and decoder.

For the language models, existing full configurations and generation metadata are
in `../../results-data/qwen-size-representation/`. The four `*-final-update.json`
files here are **extracted final records** from completed training audit logs, not
whole copies of those logs. They confirm actual exposure and one completed epoch.
The final record's `unique_source_conformers` counter is based on source records;
it is not evidence of chemically unique molecules. Table S8 therefore says records.
Precision records are the saved rank-0 audit for each run; resolved training
configuration excerpts also record gradient clipping at norm 1.0. Repeated architecture
and file listings are omitted from these excerpts; full-source hashes and the exact
selection of fields are recorded in `sources.json`.

`filter-summary.json`, `filter-per-source.csv`, `filter-normalized-summary.json`,
`filter-ref-summary.json`, and query keys document executed filtering. The ref
exclusion applies to the 1,236 matched ligands, not the alternative all-4,838 trial.
The three removal categories are assigned in order and do not double-count rows.
Rows written separately to CASF-related files were excluded from the clean training
CSV; “preserved separately” in upstream notes does not mean retained for training.
The normalization and expansion rule is verified against inspected source code;
source hashes are in `sources.json`. Two CSV receipts have only their line endings
normalized from CRLF to LF, with both original and archive hashes recorded. Normalization changes comparison keys, not
retained coordinates or molecular records. Plain and canonical-tautomer exclusions
are retained as a safety net.

Inference metadata verify sample count, temperature, top-p, seeds and decoder
binding. The inspected generation routine additionally has a 1,000-new-token limit
and filters graph mismatches, incomplete blocks and decoder failures. Historical
resolved library defaults, including top-k, were not recorded in these metadata:
we do not describe top-p=1 as unrestricted categorical sampling. Code hashes are
inspection provenance, not proof that every historical job used that exact revision.

`posebusters-config.json` is extracted from the inspected explicit configuration.
It is not an assertion of identical package versions or overrides in every run.
Historical PoseBusters membership is preserved; it was not recalculated here.

## Remaining author/reproducibility TODOs

1. **Data provenance:** original ref archive source/release/preparation and upstream
   PLINDER extraction/cleaning, including protonation, tautomers and alternate locations.
   Existing visible TODOs remain in Methods and the selection appendix.
2. **Training exports:** explain the 3,562,164-record difference between binned and
   FSQ exports. The final packaging conversions themselves report zero dropped rows;
   the discrepancy is upstream, so it must not be called tokenizer rejection without evidence.
3. **Overlap:** final residual audit against core94/ref1044/PLINDER22, including
   tokenizer development data (reused from the original validation collection).
   Record the older step-64,000 model's own training/exclusions separately. The
   executed v2 exclusion procedure is now verified; do not revert to a blanket claim
   that its similarity definition or PLINDER policy is unknown.
4. **Runtime:** fully resolved generation settings and software environments for all
   evaluated pipelines; old Qwen and external checkpoint hashes; historical RMSD
   symmetry limits, rotatable-bond option, fallback counts and failed-pair counts;
   PoseBusters versions, energy-check seeds and overrides. Current adapter defaults
   or an installed environment alone cannot resolve these.
5. **Release/citation:** source release identifiers and redistribution terms,
   permanent deposits for analysis data/code and evaluated checkpoints, final public
   CoordToken bibliographic record. A visible Data and Software Availability TODO
   distinguishes these materials from the future large conformer dataset.

Older working Markdown drafts and dated audit notes are historical. The active
LaTeX files above govern this revision. `../../METHODS_NOTES.md` records this update.

## Full deduplication account supplied 6 November 2026

`deduplication-author-description.md` preserves the supplied description; its relative
links refer to the source project, not this paper repository. The concise scientific
account is incorporated into Methods, with corpus splitting, exact-string deduplication,
filter expansion, validation scope and source-specific removal rates in the training
appendix. No additional result table or figure was introduced.

The release-status receipt, split policy, post-deduplication audit summary, validation
filter receipt and shuffle summary are archived here. Summary extractions record the
full-source hash and selected fields in `sources.json`. `qwen-validation-status.json`
records the disabled validators in all four executed v2 configurations, with source
hashes matching the existing resolved-configuration receipts. The tokenizer's input
manifest still points to the original corpus validation file; the author's description
of a separately cleaned validation export does not prove that tokenizer monitoring
used it. Its overlap remains a TODO.

**Interpretation:** core285 and the original PLINDER23 receive the normalized-match,
similarity-to-expanded-query, and match-of-similar exclusions. Ref1236 receives an
additional normalized-match exclusion only. This is one jointly filtered training
population, not separate training runs for core and ref. Some ref analogues may also
be caught by the core/PLINDER screen. “Exact” refers to normalized K1/K2 equality,
not exact SMILES or the full stereochemical key used for corpus splits.

The current Table 3 records Qwen's Hit@0.75 as 92.5532% on core and 91.6667% on ref
at 1,000 candidates, highest among generated methods in both panels. Table 2 gives
the same leading recovery rank at the ChEMBL-count target. Results now connects the
strong core recovery with its documented neighborhood exclusions. Discussion explicitly
avoids treating the core/ref difference as a causal comparison of exclusion policies:
the evaluation molecules differ, and no otherwise-matched training ablation is available.
External generators have not been shown to use these same exclusions.

The document's local label “CASF-2016 REF” does not identify the original 4,838-ligand
archive or establish its release. That source TODO remains. The supplied description
also does not reconcile the binned/FSQ export count difference or establish the older
step-64,000 checkpoint's lineage, so those TODOs remain.

## Additional September filtering records supplied 6 November 2026

See `../selection-audit/plinder/README.md` for the archived audit script, separate
JSONL writer and verified 285-core/23-PLINDER query list. Their roles agree with the
author's description, but this September workflow must not replace the October
normalized-key exclusion procedure established by the v2 input manifests. No new
checkpoint linkage has been established; the older step-64,000 TODO remains.
The author also supplied OpenReview records A0QEfcC2J9 and VYRidQUY2s for CoordToken.
Public-page retrieval was blocked by access/verification responses during this update;
the existing PDF-based citation remains unchanged pending bibliographic verification.
