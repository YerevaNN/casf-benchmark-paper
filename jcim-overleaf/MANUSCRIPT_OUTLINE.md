# JCIM manuscript outline

Agreed structure, revised 28 September 2026. This is the authors' plan, not a
journal-mandated section sequence. Follow WRITING_GUIDELINES.md for language
and reasoning, METHODS_NOTES.md for implementation/provenance, and
RESULTS_NOTES.md for numerical and figure sources.

## Governing perspective

The biologist's `docs/jcim_publication_report_2026_09_14.html` defines the
scientific perspective throughout the manuscript. The motivation is selecting
generators for conformer ensembles intended for biological applications.
The study evaluates recovery of bioactive conformations through observed bound
ligand geometries, with validity and sampling conditions specified. It does not
establish improved model training, docking, affinity prediction, or equilibrium
populations. Use corrected archived results rather than the report's historical
values. Keep excluded energy-window and training-recipe analyses outside scope.

## Central question

Which generators produce ensembles that recover bioactive conformations,
using experimentally observed protein-bound ligand geometries as structural
references, and how does performance depend on sampling budget, structural
tolerance, validity, and molecular complexity? Single-reference recovery measures whether an
ensemble reaches an observed geometry. Multi-reference coverage and precision
distinguish that reach from the concentration of samples near observed states.
Geometric diversity helps interpret these outcomes but is not assumed to be a
sufficient objective. Preserve this purpose through the transitions: sampling
and tolerance determine recovery of bound geometries; diversity is a possible
indicator of that recovery; size and flexibility identify difficult ligands;
multiple references distinguish coverage from the concentration of generated
samples near observations. Qwen is a finding within this evaluation, not the
paper's organizing objective. Use “bioactive conformations,” not “bioactive
compounds”: the task is conformational sampling of specified molecules.
Differences in energetic treatment belong to the tested
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
4. **Discussion (not drafted).** Interpret generator selection for ensembles
   intended to recover bioactive conformations: recovery from a large pool
   differs from frequent sampling near the observed set, and diversity alone
   is insufficient for selection. Discuss molecular complexity and the limits
   of precision relative to incomplete experimental references. Address
   overlap, checkpoint selection, unequal retained budgets, incomplete
   references, and the absence of downstream docking/affinity measurements.
5. **Conclusions (not drafted).** State the implications for evaluating and
   selecting bioactive conformer generators under specified sampling and
   matching conditions. Keep model comparisons within that scope, without
   extending the claims to thermodynamics, downstream utility, or universal
   superiority.
6. **Abstract and title.** Lead with evaluation of bioactive conformation
   recovery and its dependence on sampling and ensemble properties. Finalize
   after the findings and limitations are fixed.
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
