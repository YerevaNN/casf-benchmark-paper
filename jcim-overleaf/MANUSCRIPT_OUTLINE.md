# JCIM manuscript outline

Agreed structure, revised 28 September 2026. This is the authors' plan, not a
journal-mandated section sequence. Follow WRITING_GUIDELINES.md for language
and reasoning, METHODS_NOTES.md for implementation/provenance, and
RESULTS_NOTES.md for numerical and figure sources.

## Central question

Which generators produce useful ensembles of bioactive conformations, and how
does their ranking depend on sampling budget, structural tolerance, validity,
and molecular complexity? Single-reference recovery measures whether an
ensemble reaches an observed geometry. Multi-reference coverage and precision
distinguish that reach from the concentration of samples near observed states.
Geometric diversity helps interpret these outcomes but is not assumed to be a
sufficient objective. Differences in energetic treatment belong to the tested
pipelines and are not isolated as a mechanism in the current study.

Use core94 as the main cohort, ref selectively for supporting analyses, and
Qwen 1.7B FSQ step47023 as the main Qwen representative. This is a contextual
ranking, not a universal winner or a controlled model-scaling study.

## Sections

1. **Introduction (`introduction.tex`).** Begin with why molecular modeling
   needs conformer ensembles and why usefulness depends on the intended task.
   Distinguish validity, diversity, recovery, reference coverage, and sample
   concentration, and explain the limited interpretation of isolated-molecule
   energetic refinement. Place the question in prior evaluation work and
   explain why sampling budget, matching tolerance, size, and flexibility
   matter before introducing the datasets. State the contribution as a
   task-dependent comparison of generators, then explain how multiple
   references separate coverage from concentration. No separate
   literature-review section is planned. Do not claim this is the first bound
   conformer benchmark.
2. **Methods (`methods.tex`).** Datasets and identity, generator/checkpoint
   definitions, candidate tiers and validity filtering, RMSD and clustering,
   multi-reference metrics, and statistical analysis.
   Comments mark unresolved provenance and unfinished analyses.
3. **Results (`results.tex`).** Recovery at both sampling targets;
   matching-threshold dependence; geometric diversity; size and flexibility;
   separate multi-reference coverage/precision. Three compact tables and four
   figures support this sequence. Do not turn placeholders into findings.
4. **Discussion (not drafted).** Interpret task-dependent generator selection,
   best-member versus typical-sample behavior, molecular complexity, and the
   limits of observed-reference precision. Address
   overlap, checkpoint selection, unequal retained budgets, incomplete
   references, and the absence of downstream docking/affinity measurements.
5. **Conclusions (not drafted).** Answer the study questions without extending
   the claims to thermodynamics, downstream utility, or universal superiority.
6. **Abstract and title.** Finalize after the findings and limitations are fixed.
   Front/end matter in the main ACS template still require author completion.

## Figures and tables

- Table 1: core94 recovery and retained counts at both the ChEMBL-count and
  1,000-candidate targets, with the stored ChEMBL3D-PB comparison shown once.
- Supporting Table S1: core mean minimum and mean ensemble median RMSD at both targets.
- Supporting Table S2: ref recovery and retained counts at both targets.
- Table 2: selected-method size strata in core and ref, with group sizes.
- Table 3: all five external methods plus Qwen in the multi-reference panel.
- Figure 1: recovery thresholds and candidate-tier sampling endpoints.
- Figure 2: clustering radii and diversity versus recovery.
- Figure 3: size and flexibility.
- Figure 4: multi-reference coverage, precision, and contrasting examples.

The main cutoffs are 0.75 Å for recovery/coverage and 1.0 Å for clustering. Curves show threshold dependence;
they should retain crossings or failures of robustness. Fixed candidate tiers
and valid-K curves answer related but different sampling questions.

## Supporting Information and unfinished work

Reserve additional Qwen checkpoints, SFT observations, all tiers/thresholds,
PB failure details, extended cohort results, and complete per-molecule data for
SI only if included in the final story. Energy-window and random-K analyses
are outside the current manuscript scope; do not reintroduce their Methods
without associated Results. Drug results require common filtering and failure
handling.
Training-overlap, checkpoint selection, and source provenance require author
records. These gaps are marked in LaTeX comments and the notes, not filled by
assumption. The initial `supporting-information.tex` contains Tables S1 and S2; additional
SI analyses and final title/author metadata remain to be assembled.
