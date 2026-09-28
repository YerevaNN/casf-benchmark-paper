# JCIM manuscript outline

Agreed structure, revised 28 September 2026. This is the authors' plan, not a
journal-mandated section sequence. Follow WRITING_GUIDELINES.md for language
and reasoning, METHODS_NOTES.md for implementation/provenance, and
RESULTS_NOTES.md for numerical and figure sources.

## Central question

Which generators recover experimentally observed ligand conformations, and
how does their ranking depend on sampling budget, structural tolerance, and
molecular complexity? Recovery is the main outcome. Geometric diversity
helps interpret it. Multi-reference coverage and precision
address complementary aspects of ensemble usefulness.

Use core94 as the main cohort, ref selectively for supporting analyses, and
Qwen 1.7B FSQ step47023 as the main Qwen representative. This is a contextual
ranking, not a universal winner or a controlled model-scaling study.

## Sections

1. **Introduction (`introduction.tex`).** Motivate experimental-conformation
   recovery; integrate the relevant classical/learned generation and evaluation
   literature; distinguish validity, diversity, recovery, and concentration;
   state the two evaluation settings and the study's questions. No separate
   literature-review section is planned. Do not claim this is the first bound
   conformer benchmark.
2. **Methods (`methods.tex`).** Datasets and identity, generator/checkpoint
   definitions, candidate tiers and validity filtering, RMSD and clustering,
   multi-reference metrics, and statistical analysis.
   Comments mark unresolved provenance and unfinished analyses.
3. **Results (`results.tex`).** Recovery and matching-threshold dependence;
   sampling effects; geometric diversity; size and flexibility;
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

- Table 1: selected-panel recovery; core counts and minimum RMSD; ref recovery.
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
assumption. A full SI document has not yet been assembled.
