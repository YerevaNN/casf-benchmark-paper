# Methods: evidence and author follow-up

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
