# JCIM manuscript outline

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


## Ten-method PLINDER comparison, 26 October 2026

Tables 1/2/4/S4 and Figures 1/5 now include RDKit minimized and torsion minimized.
Results and Discussion distinguish coverage from precision after minimization:
primary-cutoff coverage improves, tighter-cutoff coverage declines, and minimized
torsion leads precision with its reference-seeded qualification. Methods/S8
separate the paired PLINDER optimization from separately sampled CASF pools.
The PLINDER-first framing and approved wide tables remain unchanged. Earlier
raw-only/eight-method statements are superseded by this update.

## Completed core/ref figures, 26 October 2026

Figures 3 (diversity/energy), 4 (two targets), and B1 (sampling curves) now show
CASF core and ref in rows, preserving the approved style. Ref energy and sampling
results are fully calculated using the core protocol. Results emphasizes shared
energy-spread patterns, outlier sensitivity, and the dataset-dependent ordering
for 80% recovery. Methods, captions, and Supporting Table S9 cover both sets.
PLINDER remains the lead multiple-reference evaluation; its sampling budget
analysis remains excluded. This supersedes earlier core-only scope statements.
See [current inputs and reproduction](results-data/core-ref-protocol/README.md).

## Layout and size/flexibility scope, 22 October 2026

The author selected refined wide option A for all main tables. Table 3 retains
paired targets and excludes PLINDER. Table 4 separates reference agreement
from pool checks and leaves supplied counts unranked. Figure 5 now places
PLINDER, CASF core, and CASF ref in rows, with size and flexibility in columns.
PLINDER uses COV-R; CASF uses Hit. Energy and sampling figures remain core-only.
See [Figure 5 sources](results-data/size-flexibility/README.md).

## Current scope, 21 October 2026

This supersedes the earlier PLINDER-23 display instructions below.
Tables 1/2 keep wide PLINDER/core/ref sections; Table 3 remains CASF-only.
The PLINDER panel now has 22 molecules and 1,067 cleaned references; full-pool
cluster occupancy is available. Tables 2/4/S4 share one distance source.
Figure 1 combines PLINDER diversity/proximity, COV-R and COV-P. Figure 2
compares CASF core/ref diversity/proximity. Core energy, sampling and size
follow, using unchanged measurements. PLINDER SI and review candidates are
refreshed. Results describes trends, with PLINDER leading the abstract and
Introduction as well. Qwen remains the original checkpoint pending retraining.
See [current data and reproduction](plinder22-data/README.md).

## Current scope, 19 October 2026

PLINDER's multiple PDB-derived bound conformations now lead the framing.
CASF core94 and ref remain complementary single-reference evaluations.
This supersedes the core-first scope and ref-in-SI decisions below.

- Table 1: wide PLINDER/core/ref groups, each with cluster counts at 0.5/1.0 Å
  and largest-cluster fraction at 1.0 Å. Retained counts stay in source records.
- Table 2: the same groups, with COV-R at two cutoffs and MAT-R for PLINDER,
  and Hit at two cutoffs and Best RMSD for CASF.
- Table 3: CASF core/ref only; paired candidate-target rows under each method.
  PLINDER does not have the corresponding ChEMBL-count comparison.
- Existing table method orders, stored-baseline separation, and generated-only
  best/second-best emphasis are preserved. Figures and Table 4 remain unchanged.
- Results leads with the multiple-reference question, follows the existing
  diversity/recovery/energy/sampling/size/precision sequence, and describes trends
  without repeating table values. Existing CASF figures are explicitly core-only.

The integrated values are the approved previews archived in
`results-data/three-dataset-tables/`. They retain PLINDER-23, the old Qwen checkpoint,
and current filtering differences. PLINDER-22 and retrained results are pending.
The author's reported exclusion policy is core exact/similar removal at a 0.4
fingerprint-similarity cutoff and ref exact-match removal. These changes do not
retroactively describe the archived checkpoint.

The following sections record earlier decisions; the scope above takes precedence.

## Intended scope, 6 October 2026

The benchmark is intended to guide conformer-generator selection for
structure-based drug design broadly. Pharmacophore-conditioned generation is
one potential use of the proposed resource, not its principal organizing
purpose. Keep the distinction between aligned ligand-conformation coverage
and improvements in downstream docking, screening, or model training.

## Current author sequence, 5 October 2026

The working Markdown now governs the Introduction, Methods, Results, and
Discussion. It supersedes the older recovery-first ordering below. Results
builds from geometric variation to CASF single-reference recovery and the joint
diversity/recovery plot, followed by energy dispersion and larger candidate pools,
size and flexibility, PLINDER multiple-reference coverage and precision,
and the planned conformer resource. The softer opening establishes the
purpose of this progression before introducing Table 1.

- Main Methods: scientific procedures; dataset diagram and generator summary
  are now Supporting Figure S1 and Table S8.
- Results: ChEMBL-count diversity (Table 1), recovery (Table 2), joint plot
  (Figure 1), energy dispersion (Figure 2), two-target recovery/diversity
  (Table 3 and Figure 3), size/flexibility (Figure 4), multiple references
  (Table 4 and threshold-dependent coverage/precision in Figure 5).
  Table 1 compares 0.5 and 1.0 Å cluster counts and largest-cluster occupancy
  at 1.0 Å, ordered by increasing 0.5 Å count with the stored reference last.
  Table 2 reports Hit@0.5, Hit@0.75, and Best RMSD (mean of per-entry minima),
  ordered by decreasing RMSD with the stored reference separate and last.
  Figure 1 plots 0.5 Å cluster count against Best RMSD; lower right is favorable.
  Table 3 places 1,000 candidates left and ChEMBL-count right, reporting Best
  RMSD, Hit@0.75, and 0.5 Å clusters in both groups. Order by decreasing
  1,000-candidate Best RMSD with ChEMBL3D-PB separate and last. Figure 3
  follows this order and compares Hit@0.75 and cluster count at 0.5 Å.
- Dataset plans remain in prose; the unfinished release checklist (former Table 5) is removed.
- Supporting Information: S1–S7 retain the earlier analyses and Plinder-23
  selection account; S8 describes generators and S9 reports energy summaries.
- Discussion: transferred working prose; Conclusions still pending.

Energy results use a common-hydrogen MMFF94s single-point calculation without coordinate optimization on the archived
ChEMBL-count ensembles. The cleaned Qwen evaluation and dataset release are pending.
Main figures are descriptive and show no uncertainty intervals. The author-requested
combined sampling Figure B1 follows Figure 3: panel A shows Hit@0.75 without uncertainty shading on a 50–95% axis;
both budget axes are linear. Panel B shows the 80% expected-recovery
threshold. Both panels display Qwen, LoQI, FlowR, and NExT-Mol only. Archived
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
- **0.75 Å is the primary recovery/coverage cutoff.** The author's 6 October
  revision adds Hit@0.5 to Table 2, replacing the saturated Hit@2.0 column.
  Keep 0.5/1.0 Å clustering radii distinct
  from these recovery cutoffs.
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
