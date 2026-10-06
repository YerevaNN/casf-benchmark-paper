# JCIM manuscript outline

## Intended scope, 6 October 2026

The benchmark is intended to guide conformer-generator selection for
structure-based drug design broadly. Pharmacophore-conditioned generation is
one potential use of the proposed resource, not its principal organizing
purpose. Keep the distinction between aligned ligand-conformation coverage
and improvements in downstream docking, screening, or model training.

## Current author sequence, 5 October 2026

The working Markdown now governs the Introduction, Methods, Results, and
Discussion. It supersedes the older recovery-first ordering below. Results
builds from geometric variation and the planned energy analysis to CASF
single-reference recovery, joint diversity/recovery, larger candidate pools,
size and flexibility, PLINDER multiple-reference coverage and precision,
and the planned conformer resource. The softer opening establishes the
purpose of this progression before introducing Table 2.

- Main Methods: scientific procedures; dataset diagram and generator summary
  are now Supporting Figure S1 and Table S8.
- Results: ChEMBL-count diversity (Table 1), recovery (Table 2), joint plot
  (Figure 1), two-target recovery/diversity (Table 3A/B and Figure 2),
  size/flexibility (Figure 3), multiple references (Table 4 and Figure 4).
  Table 1 compares 0.5 and 1.0 Å cluster counts and largest-cluster occupancy
  at 1.0 Å, ordered by increasing 0.5 Å count with the stored reference last.
- Dataset plans: Table 5 is explicitly an unfinished release checklist.
- Supporting Information: S1–S7 retain the earlier analyses and Plinder-23
  selection account; S8 describes generator mechanisms and training resources.
- Discussion: transferred working prose; Conclusions still pending.

Energy results, cleaned Qwen evaluation, and the dataset release are pending.
Current figures are descriptive and show no uncertainty intervals. Archived
interval calculations and older display designs remain evidence, not the
current presentation. See RESULTS_NOTES.md and `results-data/`.

## Earlier outline (28 September 2026)


Updated 28 September 2026 following the author's scope decision.

## Governing question and contribution

Which evaluated generators recover bioactive conformations under the validity
and sampling constraints defined in this study? A recovered core ligand has
at least one conformer with the requested identity that passes PoseBusters and
matches the bound geometry at **0.75 Å**. Compare the ChEMBL-count and
1,000-candidate targets, both defined before filtering, and report retained
counts. Diversity is an explanatory descriptor, not an invented pass criterion.

The main advantage relative to McNutt et al. (JCIM 2023) is breadth of evaluated
generators: four classical pipelines and six learned generators, plus the stored
ChEMBL3D-PB ensemble. Their main comparison used RDKit and DMCG to study ensemble
construction and downstream applications. Our comparison focuses on generator
selection against structural criteria; it does not demonstrate better docking
or screening. Do not describe the four classical variants as four independent
algorithm families, claim an exhaustive benchmark, or claim priority.

The biologist's HTML report remains the interpretive guide. Corrected archived
data supply the numbers. The current author instruction takes precedence over
the report's proposed threshold, training, energy, and random-K experiments.

## Scope

- **Core94 is the main CASF cohort.** Ref results remain in SI.
- **0.75 Å is the only recovery/coverage cutoff presented.** Keep the separate
  1.0 Å clustering radius clearly distinguished from recovery.
- Qwen 1.7B FSQ step47023 remains the main Qwen representative. Explain its
  selection and avoid universal-winner claims from a selected checkpoint panel.
- The 23-molecule multi-reference panel remains exploratory. Its supplied pools
  do not share a common PB filter or RMSD-failure convention.
- Do not add energy windows, random-K curves, training ablations, or downstream
  claims without new author direction and corresponding evidence.

## Narrative

1. Introduction: practical conformer-selection problem; precise relation to the
   prior JCIM study; diverse available methods; our common structural criteria;
   broader method panel; secondary multi-reference evaluation.
2. Methods: identity, datasets, generator configurations, candidate targets,
   validity checks, recovery at 0.75 Å, clustering, and uncertainty.
3. Results: core recovery and retained ensembles; diversity versus recovery;
   core chemistry-specific successes/failures; exploratory multi-reference
   coverage and concentration. Describe small differences and sparse strata
   proportionally to their uncertainty.
4. Discussion and Conclusions remain to be drafted. Connect each observation
   to generator selection under the tested constraints, preserving limitations
   on overlap, checkpoint selection, retained budget, and reference completeness.
5. Title, abstract, and author/end matter still require completion.

## Main displays

- Table 1: all ten pipelines, core recovery and mean retained counts at both
  targets, plus the single stored ChEMBL3D-PB comparison.
- Figure 1: a compact combined core recovery plot at 0.75 Å, connecting
  ChEMBL-count and 1,000-candidate endpoints for each method. No second panel;
  retained counts remain in Table 1.
- Figure 2: core diversity versus recovery for the full panel, at one clustering radius.
- Figure 3: core size/flexibility heatmaps for the full panel. Cells show exact
  recovered/total counts; the smallest groups are explicitly marked.
- Table 2: exploratory multi-reference supplied counts, COV-R, COV-P, and separate PB.
- Figure 4: multi-reference coverage/precision intervals and paired molecules,
  with coincident points counted explicitly.

## Supporting Information

- Table S1: core mean minimum and mean ensemble median RMSD.
- Table S2: ref recovery and retained counts at both targets, at 0.75 Å.
- Table S3: selected-method size-stratum values in core and ref.
- Table S4: exploratory multi-reference matching distances.

Provenance and training overlap remain unresolved author items in METHODS_NOTES.md.
Do not compile or refresh manuscript previews during ordinary source revisions.
