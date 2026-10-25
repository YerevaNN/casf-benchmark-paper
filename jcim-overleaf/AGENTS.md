# Manuscript authoring

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
