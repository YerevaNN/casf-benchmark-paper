# Methods: evidence and author follow-up

## Full corpus deduplication and panel-specific exclusions, 6 November 2026

The author-supplied full description is archived in
`method-data/training-audit/deduplication-author-description.md`. Methods now separates
full-InChIKey corpus splitting and exact-string deduplication from the later benchmark
exclusions. Core285/PLINDER23 remove normalized matches plus specified fingerprint
neighborhoods; ref1236 adds normalized matches only. All panels share the jointly
filtered training source; ref analogues may also be removed by the other screens.
Do not describe ref exact matches as literal SMILES equality or call these scaffold splits.

The training appendix covers split precedence, whole-molecule holdouts, the nablaDFT
conformation-test exception, exact-duplicate counts, the bounded similarity expansion,
training-only final filtering and separate validation treatment. The tokenizer's original
development subset still needs a residual-overlap audit. Results notes strong Qwen
core recovery under neighborhood exclusion; Discussion does not attribute core/ref
performance differences to filtering. The panels are not a controlled exclusion ablation.
The ref source, older Qwen provenance and binned-export discrepancy TODOs remain.
No benchmark scores, tables or figures changed; there are still eight SI tables and
six SI figures, included in both PDF entry points.

## Verified CoordToken and Qwen training, 6 November 2026

The active Methods now distinguishes the filtered, two-epoch CoordToken checkpoint
(37e5aa88...) from the full-corpus model in the supplied CoordToken manuscript.
The four v2 Qwen models completed one epoch, verified by their final training audit
records. The executed core/PLINDER similarity screen and ref identity exclusions
are documented, including normalized keys, ECFP4 Tanimoto > 0.4, and removal counts.
The older step-64,000 model still requires its own provenance; neither its name nor
v2 records establish its training exclusions. Residual overlap, including tokenizer
development data, remains a follow-up.

Main Methods contains the short account; `training-appendix.tex` and new Table S8
contain architecture, schedules, source counts and completed exposure.
`evaluation-protocol-appendix.tex` records checks, alignment and unresolved historical
runtime settings. Both are included through `supporting-content.tex` in the main
review PDF and standalone SI. There are now eight supporting tables and six figures;
S1–S7 and all figures retain their numbering. Results are unchanged.

Evidence and the consolidated TODO list: `method-data/training-audit/README.md`.
The CoordToken citation is now populated from the supplied PDF; its permanent
public bibliographic record remains to be confirmed. This supersedes older notes
that the v2 exclusion definition or PLINDER training policy is wholly unspecified.

## Verified filtering counts, 6 November 2026

Methods and Appendix Table S7 now account for every core/ref input using the
executed matching logs archived in `method-data/selection-audit/`.
Core: 285 minus 189 unmatched, one topology failure, and one torsional
eligibility failure leaves 94; all inputs parsed.
Ref: 4,838 minus 11 parsing/preparation failures, 3,566 unmatched, 12 topology
failures, and 13 torsional eligibility failures leaves 1,236. The frozen
identity selection removes 192 additional entries to give 1,044.
The ref archive source and upstream PLINDER cleaning remain visible TODOs.
No evaluations were rerun and no result values changed. S1–S6 remain unchanged.

## Curated Supporting Information, 5 November 2026

The 18-figure review gallery is removed from the compiled main manuscript.
The current separate SI has Tables S1–S6 and Figures S1–S6. S1/S2 figures remain
the dataset diagram and approved Qwen scatter. New S3/S4 refresh former C02/C07
for core94/ref1044; S5/S6 refresh C11/C14 for all Table 4 methods, independent
torsion pools, and the four-v2 validity comparison. Other candidates remain
historical assets only. Do not restore the review gallery during later updates.

Table S1 now records retained counts for both CASF panels/targets; S2–S4 are the
PLINDER screen, statistics, and identities; S5 is generator provenance and S6
is energy. Former S1/S3/S4 are excluded as unused or repeated summaries. Main
Methods/Results and captions use the new SI numbering. Current reproduction and
selection decisions are in `appendix-data/README.md`; use `appendix-data/build.py`
and `appendix-data/verify.py`, not the historical candidate renderer.


## Consolidated Methods, 5 November 2026

This section supersedes the older procedural descriptions below. `methods.tex`
is the current manuscript source; `method-data/methods-working.md` remains the
historical October transfer, not a current draft.

The Methods now separates panel construction, generator preparation, candidate
selection, validity, geometric metrics, recovery, energy, sampling curves,
size/flexibility, and the four Qwen configurations. Repeated generation details
were moved out of the multiple-reference metric definition. Unreported analyses,
release chronology, and proposed dataset-construction procedures were removed
from the manuscript prose. No numerical records, tables, or figures changed.

| Current result | Procedure and verified evidence |
| --- | --- |
| Core94 and ref1044 | Connectivity/stereochemistry intersection; ref selects one crystal per identity independently of model scores. `results-data/ref-panel/selection.json` and `selection.csv` retain the exact ranking, seed, and identities. |
| PLINDER22 | 1,067 cleaned references, inclusive recovery cutoffs, full supplied pools, separately pooled PB. `plinder22-data/reference-inventory.csv`, `reference-records.csv`, and nearest-distance archives. |
| Classical raw/minimized comparison | CASF uses separate raw/minimized pools and ChEMBL3D topology seeds; PLINDER pairs raw/minimized conformers. The two independent torsion runs use five graph-derived ETKDG starts with 200 accepted perturbations each. `plinder22-data/independent-torsion/` supersedes the experimental-seed descriptions below. RDKit paired minimization remains documented in `plinder22-data/minimized-runs/`. |
| Tables 1–2, Figures 1–2 | Greedy cluster counts, largest-cluster fraction, CASF hit/minimum RMSD, PLINDER COV/MAT. Benchmark `analysis/metrics.py`; current per-molecule data and cached nearest distances. |
| Figure 3 and S9 | Common-hydrogen MMFF94s single-point energies, no hydrogen or heavy-atom optimization, population SD, finite outliers retained. `energy-data/protocol.json`, `results-data/core-ref-protocol/`, and `results-data/qwen-v2/`. |
| Table 3, Figures 4/B1 | Candidate subsets precede PB; rejected candidates are not replaced. Exact uniform-subset probabilities include PB rejections in the budget. `generation/conformer_sets.py`, `generation/normalizer.py`, `analysis/common_protocol.py`, and archived sampling inputs. |
| Figure 5 | CASF-only heavy-atom and rotatable-bond strata; fixed target, all-entry recovery denominator. `results-data/size-flexibility/`. |
| Figure 6 | Connect the same molecule's MAT-R/P across nine methods; highlight the union of each method's worst molecule in either metric. `plinder22-data/molecule-tracking/`. |
| Figure 7 and Supporting Figure S2 | Four v2 Qwen runs, ref1044 and PLINDER22; paired changes, fraction strictly improved, pooled PB, and per-molecule contribution to the unweighted mean validity change. `results-data/qwen-size-representation/` archives exact run identities, representation specifications, configurations, and measurements. |

Additional checks and clarifications:

- A CASF candidate is a member of the saved pool after generator-specific
  preparation, not every generation or minimization attempt. The two targets
  select from each pipeline's pool; raw/minimized pipelines remain distinct.
- The 2,000-conformer cap in the reference loader is conditional and does not
  describe the current 1tlp evaluation: the selection manifest records 21
  conformers of its matching stereoisomer, and the current baseline retains 21.
  The irrelevant conditional cap was removed from Methods.
- Qwen binned configurations use 256 interior bins over ±12 Å and overflow
  symbols; FSQ uses 4,096 geometry codes with a molecular prefix and one code
  per atom. Config snapshots confirm a one-epoch plan, context length 4,096,
  and peak learning rate 0.0008. The two sizes share a manifest within each
  representation; training exports/token exposure differ across representations.
- The older binned 0.6B step64000 row is distinct from v2 binned step66057.
  Primary Qwen remains FSQ 1.7B step41054.
- Intended training exclusions are not an executed audit. The missing tokenizer
  citation remains commented out with the user-requested visible TODO. Training
  exclusions, full upstream reference-cleaning provenance, and alignment-fallback
  frequency remain explicit author follow-ups; they have not been filled by inference.

Validation: current table/data checks and source checks for citations, labels,
includes, artwork paths, and balanced LaTeX. No local manuscript compilation.


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


## Current ref panel: 1,044 crystals, 30 October 2026

The paper's ref panel is now the published `ref_dedup` cohort: one selected
crystal per stereochemical ligand identity, with 1,044 crystals and 1,042
ChEMBL parent compounds. Use “ref” in the manuscript, not “ref dedup” or a
sensitivity-set label. The frozen selection is independent of model outputs;
its seed, identity mapping and selection manifest are archived in
`results-data/ref-panel/`. This supersedes all 1,236-entry ref counts below.

Tables 1–3, S2/S3/S9, Figures 2–5/B1 and Supporting Figure S1 now use this
selection. Core and PLINDER data, including Table 4 and Figure 6, are unchanged.
The original Qwen 1.7B FSQ checkpoint remains in CASF comparisons. Table rankings
are recomputed from unrounded values and ChEMBL3D remains an unranked baseline.
Published common-protocol energies and exact sampling curves were imported;
no conformer generation, optimization, or alignment was rerun.

Recovery trends are similar, but the previously extreme ref LoQI mean energy
SD disappears after selection. The current Results uses Torsional Diffusion
and raw torsion as examples of high mean versus median energy spreads.
Ref energy totals are 734,206 finite of 734,207 retained conformers, with one
undefined FlowR energy and 42 empty method–entry outputs. Continuous count-target
metrics have 1,032–1,044 measurable entries; recovery uses all 1,044.

Reproduction and checks: `results-data/ref-panel/README.md`. Historical archives
and dated notes may retain the superseded denominator; they must not supply
current ref tables or figures. Figure 5 stays CASF-only.


## Additional Qwen run and molecule-level figure, 26 October 2026

Table 4 and S4 now add Qwen 0.6B CASF-dedup binned PM12 step64000 alongside
Qwen 1.7B FSQ step47023. Figure 6 follows Table 4 and connects molecule-level
MAT-R/P for the same eight methods as the approved tracking preview plus the
new Qwen run. It excludes the minimized variants as a display choice; both
remain in Table 4. The original ten-method PLINDER archive and other figures
and tables are unchanged. This supersedes blanket statements that no additional
Qwen result is included, but does not establish a leakage-free evaluation.

Table 4's renderer combines results-data/druglike_selected.csv with
plinder22-data/qwen-dedup/table-records.csv; do not overwrite the supplemental
records when refreshing the original panel. Its new values are COV-R 89.3%,
COV-P 31.8%, MAT-R 0.441 Å, MAT-P 1.113 Å, supplied N 997.8, and pooled PB 60.5%.
The new row does not alter the existing best/second-best highlights.
Figure 6 sources and the approved artwork are in plinder22-data/molecule-tracking/.
The executed training exclusion manifest, including the PLINDER policy, remains
outstanding. Different Qwen sizes and coordinate representations prevent
attributing the observed differences solely to training exclusions.


## Finalized paired PLINDER minimization, 26 October 2026

The new methods are rdkit_random_minimized_22_druglike and
torsion_minimized_22_druglike from the finalized 23 October exports. Each
minimizes the existing raw pool in its original order after rebuilding H,
using MMFF94s for up to 500 iterations, then removes H. Status 0/1 outputs
are retained, with no extra sampling, clash filter, or PB filtering of COV/MAT.
All 22 included molecules have 1,000 conformers in both runs. The manifest's
23 attempted molecules includes unsupported bortezomib, which is absent from
both generated output dictionaries and the cleaned reference panel. No UFF
fallback is used. Unconverged retained counts are recorded in the manifests.
These paired pools must not be described as separately sampled CASF-style pools.
Both torsion variants inherit the experimental seed.

The corrected raw-baseline PB evaluator renumbers atoms to the reference topology
before checking identity and stereochemistry. Its PB/energy fields are imported
only for the existing cleaned 22 identities; COV/MAT and full-pool clustering
are unchanged. See `plinder22-data/minimized-runs/` for code snapshots, hashes,
and per-molecule corrected values. No scientific scoring was rerun for this edit.

## Completed ref common-protocol analysis, 26 October 2026

Figures 3/4/B1 and Table S9 now include core/ref. Ref calculations were completed
on 22 October, using RDKit 2025.09.6, MMFF94s, rebuilt hydrogens, no optimization,
and preserved saved PB membership. Ref yielded 835,972 finite energies from
835,973 retained conformers; the one undefined FlowR energy is explicitly missing.
All finite outliers are retained. Ref energy summaries omit 50 empty outputs;
core summaries omit three. Each geometric/energy pair uses the same entries.

Sampling uses exact uniform subsets of pre-PB fixed pools for all 1,236 ref
entries, including rejections in the budget and zero-output entries in the mean.
The 63 budgets, inclusive RMSD cutoffs, and exact inversion for 80% recovery match
core. Four methods are displayed; uncertainty bands remain undisplayed.
The calculation archive includes 2,000 paired bootstrap resamples with seed
20261006. Core parity and ref full-pool endpoint checks passed upstream.
See `results-data/core-ref-protocol/README.md` for calculation sources, release,
and saved-membership reconstruction. This supersedes earlier core-only scope.

## Size/flexibility extension, 22 October 2026

Figure 5 now stratifies all three sets using the same heavy-atom/rotatable-bond
bins. CASF core/ref use archived PB-filtered fixed-target Hit@0.75; the stored
baseline is unchanged. PLINDER uses an equal-weight mean of per-molecule
COV-R@0.75 without PB filtering. No scoring or conformer generation was rerun.
Sparse groups and the empty PLINDER heavy-atom >=40 group are explicit.
See `results-data/size-flexibility/` for inputs, hashes, and reproduction.

## Cleaned reference evaluation, 21 October 2026

This supersedes earlier PLINDER counts, 200-conformer clustering, strict-cutoff
and cache/table discrepancy notes. The 20 October panel exports cover 22 molecules
and 1,067 references after bortezomib exclusion and removal of three wrong-stereo
records. Full supplied pools provide clustering/occupancy and COV/MAT; PB is
separate. All eight methods use inclusive RMSD cutoffs and have finite nearest
distances throughout. Tables 2/4/S4 and curves now agree within 1e-12.

The original Qwen checkpoint is retained. Original selection frequencies remain
archived separately from cleaned PDB/system/reference counts. Upstream cleaning
code/release lineage, protonation/tautomer and alternate-location details remain
outstanding; training-exclusion metadata is also unchanged. See
[PLINDER-22 provenance and reproduction](plinder22-data/README.md).

## Three-set table integration, 19 October 2026

Methods now describes PLINDER as the multiple-reference evaluation and includes
CASF ref in the main tables. Candidate-budget and energy analyses remain CASF-only;
the energy and existing CASF figures use core. PLINDER clustering uses up to 200
supplied conformers, whereas its coverage and matching distances use the full pool.
Table 2 reaggregates cached reference minima for all eight methods; all reference
minima are finite. The Table 4 historical summaries remain unchanged, with small
differences documented in `results-data/three-dataset-tables/cache-table-comparison.csv`.

The author reports that the revised training corpus removes core94 exact matches
and similar molecules at fingerprint similarity 0.4; ref exact matches are removed
without an additional similarity exclusion for that panel. Retraining is ongoing.
These are author-provided protocol details, not a completed manifest audit. Record
the fingerprint definition, similarity coefficient, threshold equality convention,
training-data hashes and PLINDER exclusion policy before the final evaluation.
Do not apply these claims to the old checkpoint displayed in the current tables.
PLINDER-22 results also require recomputation after the covalent-entry exclusion.

## Two PLINDER baselines added to Table 4, 12 October 2026

Added `rdkit_random_raw_druglike` and `torsion_raw_druglike` from completed
baseline COV/MAT and PB outputs. Both cover the same 23 molecular identities,
2,450 reference records, and 1,000 supplied conformers per molecule. Existing
six-method source rows remain unchanged. Sources and hashes are archived in
`results-data/druglike-baselines/` and `druglike-baselines-provenance.json`.

Torsion raw starts from each molecule's first experimental reference; it is
explicitly labeled as crystal-seeded in Table 4 and described in Methods.
RDKit removes seed coordinates before embedding. Both PB summaries report
9,000 tetrahedral-chirality failures and 1,000 double-bond-stereochemistry failures
(failure categories overlap). The torsion baseline also has 43 energy-ratio
failures. These are reported results, not corrected or silently filtered pools;
identity/reference consistency should be audited before the final evaluation.

Bold and underlining rank the best two distinct unrounded values across all
eight methods, with ties shared. Higher COV-R/P and PB pass and lower MAT-R/P
are preferred. Supplied counts are highlighted by size only, not as a quality
claim. Figure 5 remains a six-learned-method comparison. No conformers or
alignments were recalculated for this table update.

## Threshold-dependent PLINDER Figure 5, 6 October 2026

The author approved the two-panel preview for the working paper. Figure 5 now
shows molecule-averaged COV-R and COV-P against RMSD threshold, with 0.75 Å
marked. All six methods retain their established colors; no uncertainty shading
or selective emphasis was added. The earlier Figure 5 stays archived.
`plinder-threshold-data/` preserves curves, per-molecule nearest distances,
cache hashes, table comparisons, and the renderer. All 23 identities and all
reference/generated counts match table records. All nearest distances are finite.
No new generation or alignments were run. Table 4 stays unchanged. The previously
identified cache/table discrepancies remain to be reconciled before final evaluation;
they are below 0.014 percentage points, and plotted one-decimal summaries agree.

## Sampling alternatives requested for review, 6 October 2026

The current Figure B1 combines two panels and displays only Qwen, LoQI,
FlowR, and NExT-Mol. Panel A shows expected Hit@0.75 against budget without shading; its y-axis
spans 50–95%. Both budget axes are linear and have direct endpoint labels. Panel B
shows the RMSD threshold for 80% expected recovery. This replaces the two
separate B1/B2 displays without changing their calculations or main figure numbers.
Both use exact hypergeometric probabilities for uniform subsets of the saved
pre-PoseBusters pools. Rejected candidates consume budget; K is capped at N.
They are conditional finite-pool analyses, not new generation or timing tests.
The author approved this display. Sampling procedures are described in Methods
without a formula; Results reports the observations. Archived bootstrap intervals
remain available but are not displayed. No universal winner is inferred.

All 884,399 retained RMSDs were recalculated without optimization and checked
against archived per-entry minima and medians. Full-pool recovery agrees with
Table 3; source SDF hashes and candidate counts are archived in sampling-data/.
The 940 pools have 939,646 candidates before PoseBusters, including 55,247
rejections. Seventeen pools are shorter than 1,000; two have no retained output.
The old random-valid-K tables and other Qwen checkpoints were not reused.

The tolerance changes the high-budget ranking: LoQI has Hit@0.5 of 83.0%,
versus Qwen's 78.7%, while Qwen leads at 0.75 and 1.0 Å (91.5% and 97.9%).
At 1,000 candidates, 80% recovery requires 0.479 Å for LoQI and 0.506 Å for
Qwen. These are descriptive differences; band overlap is not a paired test.

## Completed common-hydrogen energy analysis, 6 October 2026

The new ChEMBL-count energy comparison supersedes the pending-energy notes
below. The archived routine used the supplied hydrogen representation and
sometimes a saved minimization energy. Inspection found no explicit hydrogens
in Qwen's retained records, while other methods generally supplied them.
The new calculation therefore removes/rebuilds hydrogens for every conformer,
evaluates single-point MMFF94s energies without optimizing any coordinates, and
ignores saved energies. This supersedes the hydrogen-relaxed protocol in commit
55ee5fa, which the author rejected in favor of raw geometries.
All heavy atoms are fixed and checked for zero displacement. The population
energy SD is calculated within each molecule and then summarized by its median
across molecules. Means and paired differences from ChEMBL3D are in Table S9.

The full protocol, RDKit version, calculation-status records, input hashes, and
per-conformer/per-molecule outputs are archived in `energy-data/`. The optimized
reference loader batches coordinate reads while preserving atom and stereo
checks. All 94 ChEMBL3D entries had zero archived PB rejections; their counts
and archived energy fingerprints were verified before rescoring. This is new
energy scoring of existing retained pools, not new conformer generation or
the pending Qwen benchmark-exclusion evaluation. The primary MMFF94s citation
is Halgren (1999), verified against the publisher's record.

## Simplified Methods and generator audit, 6 October 2026

Main Methods now describes the scientific procedure without SMARTS syntax,
function names, optimizer return codes, or validity-mask implementation.
`method-data/implementation-audit-2026-10-06.md` preserves the earlier detail,
including numerical PoseBusters tolerances and seeds, the optional reference
cap, saved-energy reuse and missing mean-energy export, and alignment fallbacks.
The corrected energy analysis remains planned, not completed.

The dataset figure and generator table move to Supporting Figure S1 and
Table S8. The table is filled from original papers/repos and local adapters;
`method-data/generator-sources.md` records the evidence and its limits.
NExT-Mol's evaluated DMT-L adapter has `use_llm=False`; ZINC/MoLlama pretraining
must not be attributed to that model. FlowR launchers specify the local
`flowr_root_v2.2_mol.ckpt`. Public `flowr_root_v2.2.ckpt` is a different,
joint-model designation; its full complex-training recipe is not inferred for
the local ligand-only checkpoint. Only published ligand-pretraining resources
are summarized as context until its exact training lineage is confirmed.


## Working Methods and appendix transfer, 5 October 2026

The current `methods.tex` follows the author's working Markdown draft, archived
with its source map in `method-data/methods-working.md`. It supersedes the older
scope description below, including a clearly marked pending energy analysis.
Plinder-23 release and identity-selection provenance are now established:
Supporting Tables S5–S7 report the screen, statistics, and all 23 identities.
Reference-coordinate preprocessing remains unresolved. Table 1 describes the
generators; Figure 1 shows the experimental panels. The initial Methods transfer retained the older Results and intervals. The
subsequent Results transfer on the same date replaces them with the working
Markdown tables and figures, which show no intervals; bootstrap provenance
below remains an archived analysis.
The Qwen tokenizer citation remains commented out with the requested visible TODO.

Both documents were compiled in an isolated directory with Tectonic and the
BibTeX fallback for local layout checks. The committed source retains Biber;
the historical preview is unchanged. The fallback emits bibliography rerun
warnings, but all cited keys and reference labels resolve in the source and
the PDFs contain no unresolved-reference markers or overflowing text. CSVs
and PDF/SVG/PNG artwork are included.

## Earlier audit (28 September 2026)


Revised 28 September 2026 against the corrected 24 September analysis.
The main panel uses one Qwen checkpoint (1.7B FSQ step47023). The main CASF figures now use only core and recovery at 0.75 Å; ref is
reported in SI. Threshold/radius sensitivity prose and the unreported ref
overlap-sensitivity procedure were removed from Methods. The diversity radius
remains 1.0 Å. Methods is now
limited to analyses used in the current Results and tables/figures. Removed
energy-window sensitivity, the unused dynamic tier, normalized entropy and
clusters per 100, drug random-K subsampling, and unused paired RMSD statistics.
Energy minimization and PoseBusters energy-ratio checks remain because they
define the evaluated generators and filtering, not a separate energy experiment.
The size export is in `docs/publication_tables_2026_09_28/`; its script
`docs/publication_size_analysis_2026_09_28.py` reads the database without mutation.
Historical checkpoint inventory below is retained for author/SI reference.
The biological framing follows `docs/bioactive_conformer_benchmark_analysis.md`;
its older counts and causal interpretations are not treated as current evidence.
This is an author working note, not part of the compiled paper. `methods.tex`
describes the implemented benchmark; it is not yet a complete account of model
training or dataset curation. Resolve the items below before submission.

## Confidence-interval provenance

The intervals are calculated secondary statistics, not dashboard-supplied
uncertainty and not repeated generation/training experiments:

- CASF paired differences: `docs/publication_analysis_2026_09_24.py`, 3,000
  bootstrap samples of ChEMBL parent compounds, retaining their complex entries.
  Exports: `docs/publication_tables_2026_09_24/paired_cluster_bootstrap.csv`.
- Drug paired differences: the same analysis script, 10,000 paired molecule
  bootstrap samples. Export: `druglike_paired.csv` in the same directory.
- Figure 4 mean-coverage bars: calculated during figure preparation by
  `manuscript/build_figures.py`, 10,000 resamples of the 23 molecule records.
  Export: `manuscript/figure-data/drug_coverage_intervals.csv`.

All are percentile intervals (2.5th and 97.5th percentiles). The scripts
initialize NumPy with seed 20260924. Marginal mean intervals and paired
method-difference intervals answer different questions. They do not account
for checkpoint selection, training randomness, or benchmark representativeness.
The user asked for the origin of the intervals; they have been retained and
their origin made explicit, not described as values read from the dashboard.

## Details still needed

1. **Dataset provenance.** Confirm the exact release and construction of the
   local `CASF16_REF` collection. Its directory name alone does not establish
   that it is the complete PDBbind 2016 refined set. For the 23-molecule set,
   recover the selection criteria, structure identifiers, source release,
   treatment of alternate locations, protonation/tautomers, stereochemistry,
   and duplicate structures. Local notes associate it with PLINDER-derived
   conformers, but the supplied pickle does not establish all these steps.
   Add the verified original dataset citation. Do not claim that these are all
   available bound conformations or independent observations.
2. **Qwen reproducibility.** Supply checkpoint hashes, training configurations,
   exact datasets/splits, tokenizer and coordinate-decoder descriptions,
   optimization settings, and inference temperature, top-p/top-k, seeds,
   retry limits and invalid-output handling. The run catalogue identifies
   checkpoints, but older descriptions of the FSQ training recipes conflict
   with the newer run labels. The draft therefore does not infer a detailed
   training history from a filename. Add the Qwen and representation references
   when that description is finalized.
3. **External inference metadata.** Archive exact checkpoint hashes, dependency
   versions and executed configurations for all five external generators.
   Inspect their preprocessing and resampling policies together. The FlowR
   scripts establish a ligand-only v2.2 model, harmonic graph inpainting,
   strict stereochemistry checking and 200 cleanup iterations by default;
   CASF permits 20 sampling rounds and druglike permits 40. Cleanup errors
   retain the unrelaxed structure. The other adapters and historical run logs
   still need a consolidated, versioned Supporting Information table.
   Defaults and current launchers should be checked against archived execution
   settings; do not describe all imported outputs as untreated samples.
4. **Druglike failure policy.** The external evaluator excludes all-nonfinite
   rows/columns from coverage denominators. The inspected 3DMolGen Qwen
   evaluator uses `nanmin`, then a Boolean threshold over all minima, making
   undefined minima misses for coverage. Both omit undefined minima for MAT.
   Thus, a shared metric name is not evidence of identical failure handling.
   The corrected audit reports nonfinite pairs for external methods. Recompute
   all methods with one policy, preserve molecule-level accounting, and revise
   the Methods when that is done. Also report raw and commonly PB-filtered
   coverage separately before claiming validity-adjusted rankings. Druglike
   random-K curves use the local finite-match implementation.
5. **CASF RMSD fallback.** `best_aligned_rmsd` falls back from symmetry-aware
   matching to identity atom indices after an exception. Its use is not counted
   in the exported results. Audit that fallback and atom correspondence before
   claiming that every score uses a verified symmetry-aware mapping. Do not
   silently remove this qualification from the Methods.
6. **Unfinished analyses.** Corrected CASF random-K curves, uniform training-set
   overlap audits and structural explanations of rescues/failures are not
   completed by this draft. Historical K curves are deliberately not described
   as corrected experiments. The core-only binned 0.6B step-29600 dynamic/count
   tiers have eight saved-count discrepancies; they are outside the shared
   core/ref panel. Confidence intervals do not account for checkpoint selection.

## Source map

Paths below are relative to the repository root unless absolute.

| Draft content | Evidence inspected |
| --- | --- |
| Cohort counts, overlap and reported comparisons | `data/mapping/casf16_{core,ref}_chembl3d_exact_intersection.csv`; `docs/corrected_findings_2026_09_24.md`; `docs/publication_tables_2026_09_24/` |
| Matching and stereoisomer resolution | `docs/data_preparation.md`; `src/casf_benchmark/chembl3d/identity.py`; ChEMBL3D loader and mapping preparation code |
| Classical sampling, MMFF94s, common PB configuration and mask projection | `src/casf_benchmark/generation/conformer_sets.py` |
| Saved-pool targets and seed 1729 | `src/casf_benchmark/generation/normalizer.py`; `scripts/run_casf_batch_item.sh` |
| Classical launcher using default seed 42 | `scripts/run_casf_ref_conformer_molecule.sbatch`; generation parser |
| Checkpoint and cohort identity | `src/casf_benchmark/config/generation_runs.yaml`; `src/casf_benchmark/config/casf_generation_families.yaml` |
| CASF scoring, conditional metrics and reference cap | `src/casf_benchmark/cli/analyze_conformer_sets.py`; `src/casf_benchmark/analysis/metrics.py` |
| Druglike PB and imported coverage results | `scripts/eval_druglike_conformers.py`; `scripts/build_druglike_covmat.py` |
| External coverage and random-K definitions | `src/casf_benchmark/analysis/druglike_covmat.py`; `src/casf_benchmark/analysis/druglike_k_efficiency.py` |
| Qwen coverage failure handling | `/mnt/weka/vtarasov/code/3DMolGen/src/molgen3D/evaluation/utils.py`, `rdkit_utils.py`, and `run_eval.py` |
| FlowR preprocessing and launch settings | `/mnt/weka/vtarasov/code/flowr_root/flowr/gen/generate_conformers_from_smiles.py`; `scripts/generate_conformers_{core,ref,druglike}.sl` in that checkout |
| Bootstrap, stratification, denominator correction and correlation | `docs/publication_analysis_2026_09_24.py` |
| Source fingerprints for the numerical report | `docs/publication_tables_2026_09_24/source_snapshot.json` |

The external source files were read in place; their full historical run lineage
has not been audited. They are not included in the Overleaf archive. Neither
the manuscript edit nor these notes changes the analysis code or results.

## Shared Qwen checkpoint inventory

The ten Qwen entries common to core and ref are:

| Run label | Representation | Parameters |
| --- | --- | --- |
| `qwen_1p7b_4e_step29600` | Binned | 1.7B |
| `qwen_1p7b_bigdata_step74000` | Binned | 1.7B |
| `qwen_1p7b_4e_from_bigdata_step26400` | Binned | 1.7B |
| `qwen_4b_4e_step20000` | Binned | 4B |
| `qwen_0p6b_fsq_6e_from_bigdata_step28320` | FSQ | 0.6B |
| `qwen_1p7b_fsq_6e_from_bigdata_step28320` | FSQ | 1.7B |
| `qwen_4b_bigdata_step110000` | Binned | 4B |
| `qwen_4b_4e_from_bigdata_step20000` | Binned | 4B |
| `qwen_0p6b_fsq_bigdata_step70534` | FSQ | 0.6B |
| `qwen_1p7b_fsq_bigdata_step47023` | FSQ | 1.7B |

The evaluated druglike panel additionally contains
`qwen_0p6b_4e_from_bigdata_step29600` and
`qwen_0p6b_bigdata_step111000`. A listed paired step-22000 run has no evaluated
rows and is not counted. Legacy core-only checkpoints are not members of the
shared panel. Transfer this inventory to Supporting Information with full
configuration identifiers before submission.

## Reference checks

Eight original method/dataset references are supplied in `acs-template.bib`.
Author spellings and bibliographic records were checked against original
papers, author repositories and DOI registration metadata. In particular,
the NExT-Mol paper is by **Liu et al.**, despite the older local catalogue's
attribution to Qiao et al.

- [CASF paper](https://doi.org/10.1021/acs.jcim.8b00545), also listed by the
  [authors' institution](https://sioc.cas.cn/sourcedb/cn/lw/202306/t20230621_6785217.html).
- [ETKDG extension](https://doi.org/10.1021/acs.jcim.0c00025), checked against the
  [authors' publication list](https://riniker.ethz.ch/publications-and-awards/journal-articles.html).
- [PoseBusters](https://doi.org/10.1039/D3SC04185A) and its
  [author repository](https://github.com/maabuu/posebusters).
- [LoQI/ChEMBL3D preprint](https://doi.org/10.26434/chemrxiv-2025-k4h7v);
  DOI registration supplied its five authors and publication date.
- [Torsional Diffusion](https://arxiv.org/abs/2206.01729).
- [MCF conference paper](https://proceedings.mlr.press/v235/wang24q.html).
- [NExT-Mol](https://arxiv.org/abs/2502.12638).
- [FLOWR.root, version 6](https://arxiv.org/abs/2510.02578v6).

The LoQI and FLOWR.root entries cite the verified preprints. Check for final
journal versions when preparing the submission bibliography. These citations
identify resources; they do not document our exact inference settings.
