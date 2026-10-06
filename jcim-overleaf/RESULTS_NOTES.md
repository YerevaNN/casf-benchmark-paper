# Results: sources and author follow-up

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

## Removed redundant displays, 6 October 2026

The author removed the pending-release checklist (former Table 5). Its unused
TeX file remains historical. Figure 5 is no longer included in Results; Table 4
retains the multi-reference results. The figure definition and artwork remain
archived. A threshold-dependent coverage/precision replacement is under discussion,
not yet included or computed for publication.

All six selected methods have cached 23-molecule RMSD matrices in the analysis
repository at `data/results/cache/druglike_rmsd_matrices/`. Reaggregating these at
0.75 Å gives differences from Table 4 of at most 0.0131 percentage points in
precision and 0.0046 percentage points in reference coverage. Resolve cache
provenance and these small discrepancies before publishing a threshold sweep;
the cache availability does not itself resolve common validity filtering.
No conformers were regenerated, no alignments rerun, and Table 4 is unchanged.

Tables 1–3 now exclude ChEMBL3D from best/second-best highlighting. Values and
row order are unchanged; tied generated values share the same emphasis.

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

## Sampling-budget revision, 6 October 2026

Table 3 now presents the 1,000-candidate target on the left and ChEMBL-count
on the right, with Best RMSD, Hit@0.75, and mean clusters at 0.5 Å for each.
Rows are sorted by decreasing unrounded 1,000-candidate Best RMSD, placing
Qwen immediately above the separate stored reference. The earlier A/B table
layout is superseded. The existing best/second-best highlighting is retained.
Figure 3 uses the same order, threshold definitions, and method palette as the
other current plots. Its connectors join measured endpoints only.

The new 0.5 Å summaries are extracted from the same pinned database, with all
recovery, RMSD, and retained-count summaries verified against casf_selected.csv.
The exported per-entry data and script are archived in results-data/ alongside
sampling-budget-provenance.json. No generation, clustering, or energy scoring
was rerun. Qwen's cluster count increases from 55.3 to 284.5; FlowR (298.1)
and random torsions (296.8) have more clusters at 1,000 candidates, while Qwen
has the lowest Best RMSD (0.307 Å) and highest Hit@0.75 (91.5%).

## Energy results and revised sequence, 6 October 2026

Results now follows geometric diversity, recovery and the joint plot, energy,
then the larger candidate budget. The energy paragraph reports observations;
it does not equate low energy SD with thermodynamic stability or absence of
implausible conformers. Figure 2 pairs mean 0.5 Å cluster counts with the median
per-molecule energy SD from the common-hydrogen recalculation. Later figures
are now 3 (budget), 4 (size/flexibility), and 5 (multiple references).
Supporting Table S9 reports mean/median SD, median per-molecule mean energy,
and median paired mean-energy differences from the same molecule's ChEMBL3D
ensemble. Every method retains its original measurable-entry set.

Qwen and FlowR have median SDs of 17.6 and 11.6 kcal/mol, respectively; random
torsion sampling and Torsional Diffusion have 40.5 and 26.5. LoQI has 3.7,
minimized baselines have 2.7–2.8, and ChEMBL3D 2.7. Torsional Diffusion's mean
SD is 1,628.4, and Qwen's is 49.0, showing the influence of extreme entries.
The paired mean-energy difference from ChEMBL3D is +20.7 for Qwen and +0.8 for
FlowR. These single-point calculations use rebuilt, unoptimized hydrogens and
supersede both the original mixed-hydrogen energies and the hydrogen-relaxed
calculation in commit 55ee5fa. No atom is optimized during energy evaluation.
See METHODS_NOTES.md and energy-data/ for the audit.

## Tighter recovery cutoff, 6 October 2026

The author replaced Table 2's nearly saturated Hit@2.0 column with Hit@0.5.
The current metric order is Hit@0.5, Hit@0.75, Best RMSD; generator row order,
reference separation, and highlighting conventions are unchanged. Both hit
rates were verified against per-entry RMSD minima using all 94 entries.
LoQI leads Hit@0.5 (67.0%), followed by ChEMBL3D (63.8%), while Qwen reaches
56.4%. The Results now describes that rank change and explicitly notes that
the other generators with the highest cluster counts, FlowR and random
torsion sampling, have higher Best RMSD. This does not assert that every
diverse generator performs poorly or that all Qwen samples are close.
Figure 1 adds dashed horizontal and vertical guides at the stored ChEMBL3D-PB
values (21.7128 mean clusters, 0.52538 Å Best RMSD); these are reference guides,
not an equality line between quantities with different units or a fitted trend.

## Recovery percentages and diversity--RMSD plot, 6 October 2026

Table 2 now contains three metrics: Hit@0.75, Hit@2.0, and Best RMSD.
The next larger threshold actually archived in the selected CASF snapshot is
2.0 Å, not 1.0 or 1.5 Å. Both percentages were verified from the per-entry
minimum RMSDs against all 94 entries, with undefined outputs counted as misses.
Best RMSD is the existing mean of those minima, verified without changing its
denominator (92 NExT-Mol, 93 MCF, 94 others). It is not the mean distance of
all generated conformers or a single panel-wide minimum. Generators run from
highest to lowest unrounded RMSD, then ChEMBL3D-PB follows a dashed rule.

Figure 1 now shows 0.5 Å cluster count on x and Best RMSD on y. The favorable
direction is lower right. All points use unchanged archived measurements;
axes show the observed region and no uncertainty intervals are added. Qwen
combines the largest mean cluster count with the lowest Best RMSD, but the
0.019 Å difference from ChEMBL3D is modest. The Results describes that joint
position without asserting statistical outlier status or that every generated
conformer is close. Archived points and the revised renderer are in results-data/.
The initial Qwen checkpoint and pending benchmark-excluded evaluation remain
explicit. The source database digest still matches the corrected snapshot.

## Cluster counts and occupancy, 6 October 2026

Table 1 orders generators by increasing unrounded mean cluster count at 0.5 Å,
with Qwen last among generators and ChEMBL3D-PB in a separate final row below
a dashed rule. Both 0.5 and 1.0 Å counts remain. Bold and underlining identify
the best and second-best distinct values in each column, including tied values
among generated methods only; the stored reference is unranked. Retained-count highlighting concerns yield only.

The added column is the mean largest-cluster fraction at 1.0 Å, expressed as
a percentage. Smaller values indicate less concentration in one cluster.
This directly interpretable occupancy measure complements cluster count;
it does not imply identical rankings or biological relevance. The database
digest matches the corrected snapshot. Read-only extraction confirms that
all three cluster statistics share identical measurable-entry sets within
each method (92 NExT-Mol, 93 MCF, 94 others) and reproduces all prior counts.
The query and aggregation conventions are archived with the updated CSV.

## Display relocation, 6 October 2026

Methods' dataset diagram and generator summary moved to Supporting Figure S1
and Table S8. Main Results tables are now 1–5 and figures 1–4. Stable labels,
CSV values, and artwork filenames are unchanged. The display numbers in the
historical audit sections below refer to earlier manuscript layouts.


## Two clustering radii in Table 2, 5 October 2026

Table 2 now reports mean clusters at 0.5 and 1.0 Å, sorted by the unrounded
0.5 Å mean. Both columns use identical measurable entries per method: 92 for
NExT-Mol, 93 for MCF, and 94 for the others. The read-only database digest
matches the corrected September snapshot, and every recomputed 1.0 Å mean
matches the previous table. The CSV and query provenance are in `results-data/`.
The text describes persistent differences in geometric diversity despite rank
changes, not statistical significance or a ranking of overall conformer quality.
Other figures and Table 4 retain the 1.0 Å radius.

## Current working Results transfer, 5 October 2026

The author's Markdown Results and Discussion now supply the LaTeX body. The
opening was softened to introduce the purpose of conformational variation
before the first comparison. The sequence is diversity at the ChEMBL-count
target, pending energy analysis, recovery at that target, diversity versus
recovery, the 1,000-candidate comparison, molecular difficulty, multiple
references, and the planned dataset release. Values were transferred without
recalculation. `results-data/results-working.md` preserves the edited source.

Tables 2–6 and Figures 2–5 follow that sequence; Table 4 keeps its A/B panels.
Methods still holds Table 1 and Figure 1. Discussion is in `discussion.tex`.
Table 6 is explicitly a pending-release checklist, not a measured result.
Author instructions remain source comments, while the planned status of
energy, training exclusions, and dataset release remains visible in the prose.

The current figures are the uniform working-draft figures, including the
ChEMBL-count diversity/recovery plot and descriptive multi-reference plot
without intervals. The older interval prose and figure definitions have been
replaced; the archived bootstrap calculations below remain available for later
editing. Methods no longer claims that the current displays show intervals.
The current source CSVs, artwork hashes, and unchanged original figure script
are archived in `results-data/`. Historical assets remain separate.

Local render check: isolated Tectonic compilation with the BibTeX fallback;
the committed Overleaf source retains Biber. Tables, figures, and Discussion
were visually inspected; the historical preview was left unchanged.

## Previous Results audit (28 September 2026)


Revised 28 September 2026. The manuscript now follows the agreed recovery-first
story: core94 anchors the comparison; Qwen 1.7B FSQ pretraining step47023 is
the sole main-text Qwen; ref is reported in SI; the 23-molecule
multi-reference drug panel is a separate experiment. Additional Qwen variants
and SFT comparisons are reserved for possible Supporting Information.

## Interpretation and numerical sources

The biologist's `docs/jcim_publication_report_2026_09_14.html` supplies the
interpretive perspective, especially its distinction between recovery,
sampling efficiency, diversity, and multi-reference precision. Its older
numbers do not supersede the corrected exports.

The source database fingerprint was checked against the 24 September snapshot
before the size export. No generation, energy evaluation, or CASF random-K
analysis was rerun in this manuscript revision.

| Content | Source relative to repository root |
| --- | --- |
| Recovery at 0.75 Å, retained counts, minimum/median RMSD, diversity | `docs/publication_tables_2026_09_24/all_casf_summary.csv` |
| Paired intervals and cohort-overlap sensitivity | `docs/publication_tables_2026_09_24/paired_cluster_bootstrap.csv` |
| Smaller sampling targets | `docs/publication_tables_2026_09_24/tiers.csv` |
| Size strata and Table S3 | `docs/publication_tables_2026_09_28/size_strata.csv`, reproduced by `docs/publication_size_analysis_2026_09_28.py` |
| Flexibility | `docs/publication_tables_2026_09_24/flexibility.csv` |
| Table 2 coverage and Table S4 MAT | Molecule means of `docs/publication_tables_2026_09_24/druglike_per_molecule.csv` |
| Table 2 PB | Pooled pass fraction in `docs/publication_tables_2026_09_24/druglike_summary.csv` |
| Drug uncertainty and examples | `docs/publication_tables_2026_09_24/druglike_paired.csv` and `druglike_cases.csv` |

Table 1 contains the full selected generator panel plus ChEMBL3D-PB, with
core recovery and retained counts at both the ChEMBL-count and fixed targets.
Supporting Table S1 contains core RMSD summaries for both tiers, and Table S2
contains ref recovery and retained counts. The exact 42 selected source rows
are archived in `figure-data/recovery-table-tiers.csv`; recovery and counts
agree with `casf_summary.csv` for every overlapping method/tier.
The main Results introduce both tiers at 0.75 Å before diversity and core
molecular complexity. The selected-method size table is now Table S3.
Table 2 reports supplied counts, coverage, precision, and separately assessed
PB for all five external methods and Qwen; matching distances are Table S4.

Historical threshold result (removed from the current single-cutoff manuscript): on core, LoQI exceeds Qwen at 0.5 Å
(83.0% versus 78.7%) despite Qwen's higher recovery at 0.75 Å. Therefore the
text does not claim threshold-invariant superiority. No correlation from the
old 19-generator panel is presented as a result for the reduced main panel.

## Audit of the opening recovery counts

The five opening counts were independently recounted from the public-release
SQLite `per_ligand_long` table for `ligand_set=core`, the fixed candidate tier,
and the exact selected methods. For each of 94 entries, `casf_best_rmsd <= 0.75`
was checked against `casf_hit_0p75`; all flags agreed. Counts were 86 for
Qwen, 84 for Torsional Diffusion, 81 for LoQI, 76 for RDKit raw, and 74 for
ChEMBL3D-PB. They are direct ligand counts, not bootstrap estimates or counts
of generated conformers. The source hash still matches the public release.

`figure-data/core-recovery-audit.csv` preserves each ligand's minimum RMSD,
hit flag and retained count; `core-recovery-counts.csv` contains the five
summary rows. The revised Results defines recovery in its opening paragraph.
Checkpoint steps are retained in Methods/provenance, not Results/captions.
Prose was revised throughout Results to connect observations and their
interpretation without changing the numerical findings.

## Figures recreated from the public dashboard

On 28 September, GitHub's latest public dashboard release was verified as
`dashboard-data-stereo-identity-v2` (18 September). Both database digests match
the local files exactly. Figure sources and explicit denominator/checkpoint
choices are archived in `figure-data/README.md` and `provenance.json`.

Four main figures are placed by commands in `figures.tex`:

1. `core-recovery-budget.pdf`: one compact core recovery plot at 0.75 Å,
   connecting both candidate-target endpoints for each generator. The stored
   ChEMBL3D-PB ensemble is a dashed reference; retained counts remain in Table 1.
2. `diversity-recovery.pdf`: the full main panel at the 1.0 Å clustering radius.
3. `size-flexibility.pdf`: core-only heatmaps with recovered/total counts for
   all ten pipelines and ChEMBL3D-PB; sparse groups are marked.
4. `multireference-coverage.pdf`: supplied-pool recall/precision and paired
   molecules, with the 11 coincident complete-coverage molecules annotated.

The additional core figure sources are `core-comparison.csv` and
`core-strata.csv`. Their source paths and digests are in
`core-figure-sources.json`; all overlapping archived strata agree. The full
panel reveals that MCF recovers 8/11 core ligands with 7–8 rotors versus Qwen's
6/11; this descriptive result now appears in the text. No new generation or
hypothesis test was performed. The original bootstrap draw order is preserved.

PDFs, SVGs, and 450-dpi PNGs are in `figures/`. The public release's unchanged
14 September energy sidecar was recreated as `energy-window-historical.pdf`,
clearly labeled historical and excluded from the paper. No current energy
robustness result is inferred from it. Additional Qwen labels in the release
were not substituted for the exact step47023 checkpoint.

The figure script also archives the dashboard's displayed summaries. MCF and
NExT-Mol core Hit@0.75 differs from the displayed conditional percentage
(81/93 = 87.1%) because the paper counts the missing entry as a miss
(81/94 = 86.2%). All manuscript figures retain the paper's defined denominator.

## Remaining scientific decisions

- Energy windows, dynamic tiers, normalized entropy, and drug random-K
  analysis are outside the current story; their Methods and energy-related
  Introduction/Results prose were removed. Historical artifacts remain
  archived but are not manuscript experiments.
- Confidence intervals were calculated from the underlying records: paired
  statistics in the corrected 24 September analysis, and marginal Figure 4
  intervals during figure preparation. See METHODS_NOTES.md for provenance.
- Harmonize drug PB filtering and RMSD-failure handling; update prose, tables,
  and figures together. Current values are descriptive supplied-pool results.
- Document representative-checkpoint selection, training/inference metadata,
  training-set overlap, and dataset provenance (see METHODS_NOTES.md).
- Keep sparse size strata descriptive. Size is not an independent causal
  explanation after stratifying on it alone.
- Additional Qwen/SFT and dynamic-tier analyses remain outside the current
  manuscript; no nonexistent SI table numbers are cited.

The Introduction cites prior experimental-conformation benchmarking rather
than claiming priority. The added McNutt et al. reference was verified against
https://pubs.acs.org/doi/10.1021/acs.jcim.3c01245 on 28 September 2026.

## Editorial revision

The prose now states findings directly, with headings describing the scientific
comparisons. Results retains selected quantitative contrasts and uncertainty
rather than repeating table rows. All numerical table entries remain unchanged.
Methods retains the settings and denominators required for reproducibility.
The positive interval for Qwen's precision deficit is the sign reversal of
the existing Qwen-minus-LoQI interval, not a new statistical calculation.

## Biological framing and transitions

The Introduction and Results transitions were checked against sections 2, 5,
7, 12, and 14 of `docs/jcim_publication_report_2026_09_14.html`. The guiding
question is generator selection for recovery of bioactive conformations,
represented by observed bound ligand geometries. Sampling and matching
tolerance qualify that recovery; diversity is an incomplete indicator; size
and flexibility identify difficult ligands; the drug panel distinguishes
coverage of multiple observations from concentration near them. The report's
older numerical values, training claims, and proposed additional experiments
were not imported into the manuscript. Drug reference provenance and common
filtering remain unresolved as documented above.

The user explicitly confirmed that this perspective governs the whole paper,
including the Introduction and Results. It must also guide the future
Discussion, Conclusions, abstract, and title. The report is interpretive
guidance; corrected exports remain authoritative for the numerical evidence.

## Presentation decision after reading the prior JCIM paper

McNutt et al., DOI 10.1021/acs.jcim.3c01245, full text:
https://pmc.ncbi.nlm.nih.gov/articles/PMC10647020/
Their main RDKit/DMCG comparison studies ensemble construction and downstream
performance. Our Introduction now states the broader method panel and our
structural selection criteria explicitly. Presentation follows a practical
question, defined evaluation conditions, then task-specific observations;
it does not import their downstream claims, numerical thresholds, or conclusions.
The main paper uses only core recovery at 0.75 Å; ref and the duplicate size
table are in SI. Figure PDFs were regenerated; manuscript previews were not compiled.

## Additional figure review, 28 September 2026

The author requested as many useful figure options from the HTML analyses as
possible before selecting the final set. `figure-candidates.tex` therefore adds
18 review-only figures at the end of the current Overleaf draft; these are not
18 newly endorsed main-text claims. `FIGURE_CANDIDATES.md` maps each option to
the report and numerical sources. The original four main figures are unchanged.

C01–C10 use the corrected core94 snapshot, the ten selected generation pipelines,
and the stored reference pool where relevant. C11–C18 retain the exploratory
23-molecule drug comparison and explicitly distinguish recall, precision, and
validity. No energy-window, random-K, reference-cohort, training-recipe, or
additional recovery-threshold narrative is restored. All continuous-distance
and diversity summaries preserve their measurable-ensemble denominators.
New paired confidence intervals use 10,000 resamples with seed 20260928 and
are labeled exploratory, without multiplicity or checkpoint-selection adjustment.
They are generated from archived per-entry or per-molecule rows; full provenance
and source hashes are in `figure-data/candidates/provenance.json`.
