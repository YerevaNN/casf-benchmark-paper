# Manuscript authoring

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

## Appendix placement in the review manuscript

The curated appendix is included at the end of the main PDF again. Both
`acs-template.tex` and `supporting-information.tex` use `supporting-content.tex`,
so the eight supporting tables (S1–S8) and six supporting figures (S1–S6) stay
identical in either document. Keep this appendix visible in the main review
PDF; only the obsolete 18-candidate gallery is excluded. The separate SI entry
point remains available for eventual submission.


## Verified filtering account, 6 November 2026

Methods and Appendix Table S7 now give exact core/ref exclusion counts from
archived matching logs in `method-data/selection-audit/`. Existing S1–S6 numbers
are unchanged. The original ref archive source and upstream PLINDER cleaning
remain visible TODOs; do not infer either from local directory names.

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


## Independent PLINDER torsion baselines

The current raw/minimized PLINDER torsion rows use
`torsion_independent_raw_druglike` and `torsion_independent_minimized_druglike`.
Both start from five ETKDGv3 conformers built from isomeric SMILES, without
experimental coordinates. The reference-seed dagger and qualification are removed.
Tables 1/2/4/S4 and Figures 1/6 use the new results. LoQI now leads PLINDER precision;
raw torsion slightly exceeds FlowR in 0.5 Å clusters. CASF results are unchanged.
Exact sources, seed audit, and reproduction are in `plinder22-data/independent-torsion/`.
This supersedes earlier instructions to mark these two rows as reference-seeded.

## Approved Qwen size and representation comparison

Figure 7 compares the four v2 configurations (0.6B/1.7B, binned/FSQ), with size
effects above representation effects. Supporting Figure S2 contains the approved
simplified molecule-level scatter, CASF zooms, and only Actinonin labeled.
Both use ref1044 and PLINDER22. The v2 binned 0.6B step66057 checkpoint in these
figures is distinct from the older CASF-dedup step64000 row in Table 4/Figure 6.
Do not replace existing table rows with these additional runs. Training exposure
differs between representations; do not describe this as an isolated causal
encoding effect. Sources and reproduction: `results-data/qwen-size-representation/`.

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


## Figure 5 scope, 26 October 2026

Figure 5 now shows CASF core and ref only, in two rows with size and flexibility
columns. The author removed PLINDER from this display; do not restore it when
refreshing figures. PLINDER measurements remain archived, and its other figures
and tables retain their approved scope. This supersedes earlier three-dataset
Figure 5 instructions below.

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


## Figure 1 panel layout, 26 October 2026

Figure 1 places A across the top and B/C at bottom left/right. A retains all ten
methods. B/C show only RDKit raw, MCF, FlowR, Torsional Diffusion, NExT-Mol, LoQI,
and Qwen, in that legend order. The existing palette is preserved; RDKit raw
has a dashed curve. This supersedes the earlier all-ten-method curve scope.
Tables and Figure 5 retain the minimized and reference-seeded torsion results.

## Finalized PLINDER minimized baselines, 26 October 2026

PLINDER now contains ten methods, adding `rdkit_random_minimized_22_druglike`
and `torsion_minimized_22_druglike`. Tables 1/2/4/S4 and Figures 1/5 include both.
Use the existing wide layouts and palette; raw classical curves are dashed,
minimized curves dotted. Mark both PLINDER torsion variants as reference-seeded.
These runs minimize the existing raw conformers, unlike the separately sampled
CASF pools. The paper archive records 22 x 1,000 conformers per run, MMFF94s,
rebuilt hydrogens, 500 iterations, and retention of unconverged outputs.
Raw baseline validity/energy records use corrected atom-order checks restricted
to the cleaned 22 identities; their COV/MAT/clustering remain unchanged.
Minimized torsion now leads PLINDER COV-P and MAT-P at the primary cutoff;
LoQI leads precision among learned methods. The six-learned-method review
candidates keep their stated subset. See `plinder22-data/minimized-runs/README.md`.
This supersedes earlier eight-method and raw-only PLINDER statements below.

## Completed core/ref energy and sampling, 26 October 2026

Figures 3, 4, and B1 now contain CASF core/ref rows using the approved style.
Their stable artwork names are `figure-energy-dispersion`,
`figure-3-sampling-budget`, and `figure-budget-combined`. Rebuild them with
`results-data/core-ref-protocol/build_figures.py`. Supporting Table S9 has two
stacked dataset sections, rebuilt by `build_energy_table.py` in that directory.
The ref energies are completed common-hydrogen MMFF94s calculations, not archived
energies or placeholders. Sampling curves use the same exact pre-PB protocol as
core. Preserve finite outliers and record the one undefined ref FlowR energy.
Results describes trends and the different core/ref ordering at 80% recovery.
Main tables, PLINDER, and original Qwen identity remain as previously approved.
This supersedes all earlier core-only energy/sampling instructions below.

## Approved layout refinement, 22 October 2026

The author selected option A for all four main tables: full-width typeset tables,
short headers, light alternating shading, and no resizebox shrinkage. Tables 1/2
retain PLINDER/core/ref groups. Table 3 retains paired targets per method and
CASF core/ref only. Table 4 groups reference agreement separately from pool
checks; supplied count is descriptive and no longer ranked. All values and
method orders are unchanged.

Figure 5 now contains PLINDER/core/ref rows and size/flexibility columns.
PLINDER reports mean COV-R@0.75 from full supplied pools; CASF reports Hit@0.75
from PB-filtered fixed-target pools, with the stored baseline unchanged. Sparse
and empty groups remain explicit. Its renderer and source lineage are in
`results-data/size-flexibility/`; the stable artwork filename is
`figures/figure-4-size-flexibility.*`. Energy and sampling figures remain core-only.
This section supersedes conflicting earlier layout and size-figure instructions.

## Cleaned PLINDER-22 update, 21 October 2026

This section supersedes earlier PLINDER counts and display instructions below.
All active PLINDER displays now use 22 molecules and 1,067 cleaned references
from the 20 October release. Bortezomib and three wrong-stereo reference records
are excluded. Clustering and occupancy use full supplied pools, not 200-conformer
subsets. Coverage uses RMSD <= threshold for every method; all nearest distances
are finite. PB is measured separately and is not applied to coverage/clustering.
The Qwen identity remains 1.7B FSQ step47023; no retrained outputs are substituted.

`plinder22-data/` archives per-molecule results, nearest distances, reference
records, input hashes, and reproduction scripts. Tables 1/2 retain PLINDER first,
then CASF core/ref; PLINDER largest-cluster fractions are now populated. Table 3
remains CASF-only. Tables 2/4/S4 and all curves share one consistent distance source.
S6/S7 and Supporting Figure S1 describe the cleaned collection. S5 retains the
original identity-screen criteria. Historical files named plinder23 are retained
where needed for stable includes/labels; the current companion identity CSV is
`method-data/plinder22-molecules.csv`.

Figure 1 now combines PLINDER diversity versus MAT-R, COV-R, and COV-P for eight
methods. Figure 2 compares CASF core/ref diversity versus Best RMSD. The remaining
energy, budget, and size figures are core-only, using unchanged measurements.
The old standalone PLINDER threshold figure is no longer included. PLINDER review
candidates C11-C18 are refreshed; C18 retains three examples after bortezomib removal.
CASF artwork and Table 3 are unchanged. Preserve the existing colors and wide tables.

Rebuild tables with `results-data/three-dataset-tables/build_tables.py` and
`results-data/build_druglike_table.py`; current plots/inventory with the scripts in
`plinder22-data/` and `results-data/three-dataset-tables/plot_casf.py`.
Do not compile a local manuscript preview. Validate values and source references;
the user reviews the compiled paper in Overleaf.

Outstanding upstream metadata: original cleaning/extraction code, full reference
release lineage, protonation/tautomer and alternate-location handling. Existing
Qwen training-exclusion and retraining follow-ups remain. The refreshed caches
resolve the old table/curve inconsistency, not those separate metadata questions.

This is the canonical manuscript directory in `YerevaNN/casf-benchmark-paper`.
Follow the repository-root AGENTS.md for the authorized commit/push workflow.

Before drafting or revising this paper, read:

- [WRITING_GUIDELINES.md](WRITING_GUIDELINES.md) for language, logical coherence,
  and handling of evidence.
- [MANUSCRIPT_OUTLINE.md](MANUSCRIPT_OUTLINE.md) for the scientific question,
  section structure, literature placement, and unfinished analyses.

Use these as authoring guidance, not as evidence that a result has been measured
or that a journal requirement has been verified. Follow the user's current
instructions when they change the scope or structure.

The user has explicitly designated the biologist's report as the governing
scientific perspective for the whole paper, including Introduction and Results.
Before changing the framing or transitions, consult
`/mnt/weka/mbedrosian/code/casf-benchmark/docs/jcim_publication_report_2026_09_14.html`
(especially sections 2, 5, 7, 12, and 14). The central contribution is evaluation
and selection of generators for recovery of bioactive conformations, using
observed bound ligand geometries as structural references. Sampling budget
qualifies recovery; diversity is an incomplete indicator; size and flexibility
identify difficult ligands; multiple references distinguish coverage from
concentration near observations. Qwen performance is a result within that
evaluation, not the organizing purpose. The corrected archived data supply
the numerical evidence: do not restore the HTML's older values, treat its
proposed experiments as completed, or reintroduce analyses the user excluded.
Preserve this framing when drafting Discussion, Conclusions, abstract, and title.


## Historical scope before the cleaned release

The author's scope decision (19 October 2026) leads with PLINDER's multiple
experimental references, with CASF core94 and ref as complementary evaluations.
Tables 1–2 have three wide column groups in that order. Table 1 repeats cluster
counts at 0.5/1.0 Å and largest-cluster fraction at 1.0 Å; unavailable PLINDER
occupancy remains pending. Table 2 uses PLINDER COV-R@0.5, COV-R@0.75 and MAT-R,
and CASF Hit@0.5, Hit@0.75 and Best RMSD. Retain the existing method order.
Table 3 excludes PLINDER: core/ref are column groups, with paired 1,000 and
ChEMBL-count rows per method and the stored baseline separate and last.
Its order remains decreasing core Best RMSD at 1,000 candidates.
Rebuild Tables 1–3 with results-data/three-dataset-tables/build_tables.py.
Results should describe trends without repeating rows of numerical results.

All displayed results still use the original Qwen checkpoint and PLINDER-23.
Do not relabel them as retrained or PLINDER-22 results. The author reports core
exclusions at fingerprint similarity 0.4 and ref exact-match exclusions; the
fingerprint definition, threshold boundary and PLINDER training-overlap policy
still require records. PLINDER currently uses supplied pools and up-to-200
subsets for clustering; CASF uses PB-filtered ensembles. Preserve these distinctions.

Figures are unchanged by the table revision: Figures 1–4 and B1 describe core,
and Figure 5 describes PLINDER. Figure 1 compares 0.5 Å clusters with Best RMSD;
Figure 2 uses the common-hydrogen MMFF94s energy protocol in energy-data/.
Figure 3 retains the Table 3 method order. Clustering radii remain distinct
from recovery cutoffs. See MANUSCRIPT_OUTLINE.md for current placement.

Read METHODS_NOTES.md before revising Methods: it records evidence and unresolved
metadata. methods.tex contains the benchmark Methods draft and is included by
acs-template.tex. Read RESULTS_NOTES.md before revising results.tex: the HTML
report supplies the interpretation, while corrected exports supply the numbers.
The dataset diagram and generator table are Supporting Figure S1 and Table S8.
Main Methods describes scientific procedures without code-level syntax.
introduction.tex is the Introduction draft. Tables are in tables/ and numbered
figure definitions are in figures.tex, with artwork in figures/.
The current figure artwork and source records are identified in README.md and
results-data/. The older build_figures.py and figure-data/ remain historical. Methods must describe only the analyses used in Results, tables, or figures.
Energy-window experiments remain outside scope. The author has combined the sampling analyses
after Figure 3 as Figure B1, showing Qwen, LoQI, FlowR, and NExT-Mol only.
Panel A uses Hit@0.75 on a 50–95% axis without shading; both budget axes are
linear with direct endpoint labels. Panel B gives the RMSD threshold for 80% expected
recovery. Keep Results observational and explain sampling in Methods without a formula.
They use exact uniform-subset recovery with rejected candidates retained in the
budget; see sampling-data/ and sampling-figure-options.tex. Historical valid-only
random-K analyses are not their data source. Front and end matter still contain
upstream examples; Discussion is a working draft and Conclusions remain
to be written. Replace
sample names and claims only with verified information or explicit placeholders.
For manuscript-only edits, check citations and cross-references in the source.
Do not compile local previews: the user reviews changes in Overleaf.
Python pipeline tests are not needed.

For Tables 1–3, best/second-best highlighting ranks generated methods only.
ChEMBL3D is an unranked baseline, always shown without bold or underlining.

The dataset-release checklist (former Table 5) remains removed. Figure 5 now
shows COV-R and COV-P versus RMSD threshold for the six PLINDER generators,
with 0.75 Å marked. Preserve Table 4 and the earlier figure assets as history.
The curves use cached distances and are archived in plinder-threshold-data/.
Small cache/table differences remain an explicit author follow-up before final evaluation.

Table 4 now includes RDKit raw and crystal-seeded torsion raw alongside the six
learned methods. Highlight the best two distinct unrounded values, sharing ties:
higher COV-R/P and PB pass, lower MAT-R/P. Mean supplied count is ranked by
size only. Keep the reference-seeded torsion label and Methods qualification.
Figure 5 remains the six-learned-method threshold comparison.
