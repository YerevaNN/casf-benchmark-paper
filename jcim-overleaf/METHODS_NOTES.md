# Methods: evidence and author follow-up

## Sampling alternatives requested for review, 6 October 2026

Figures B1/B2 follow the main budget figure without changing main figure numbers.
B1 shows expected recovery against candidate budget at 0.5, 0.75, and 1.0 Å,
with pointwise 95% intervals from 2,000 paired entry-bootstrap resamples. B2
shows the smallest observed RMSD threshold reaching 80% expected recovery.
Both use exact hypergeometric probabilities for uniform subsets of the saved
pre-PoseBusters pools. Rejected candidates consume budget; K is capped at N.
They are conditional finite-pool analyses, not new generation or timing tests.
The author will select the final figure later; no universal winner is inferred.

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
