# JCIM manuscript

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


## Independent torsion update

The PLINDER raw/minimized torsion rows now use independent ETKDG starting
conformations. Tables 1/2/4/S4 and Figures 1/6 have refreshed values and rankings;
the reference-seed marker and explanation are removed. Methods describes the
five-start generation procedure, and Results/Discussion reflect the new ranking.
CASF and all other methods are unchanged. See the
[source archive and checks](plinder22-data/independent-torsion/README.md).

## Qwen size and representation comparison

Figure 7 adds the approved two-row comparison: model size above and coordinate
representation below, using CASF ref matching, PLINDER matching, and PLINDER
validity. Supporting Figure S2 contains the simplified molecule-level scatter
with CASF zooms and Actinonin labeled. It is included in the separate
`supporting-information.tex` document. Results explains both; Methods distinguishes
the four v2 checkpoints from the older 0.6B CASF-dedup row in Table 4.
The existing tables and Figures 1–6 are unchanged. Inputs, exact identities,
training configurations, source hashes, and reproduction scripts are archived in
[results-data/qwen-size-representation/](results-data/qwen-size-representation/README.md).

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


## Figure 1 panel layout, 26 October 2026

Figure 1 places A across the top and B/C at bottom left/right. A retains all ten
methods. B/C show only RDKit raw, MCF, FlowR, Torsional Diffusion, NExT-Mol, LoQI,
and Qwen, in that legend order. The existing palette is preserved; RDKit raw
has a dashed curve. This supersedes the earlier all-ten-method curve scope.
Tables and Figure 5 retain the minimized and reference-seeded torsion results.

## PLINDER minimized runs, 26 October 2026

Tables 1, 2, 4 and S4 and Figures 1 and 5 include the finalized RDKit minimized
and torsion minimized PLINDER-22 runs. Methods, Results, Discussion, captions,
and S8 now distinguish paired minimization of the PLINDER raw pools from the
separate CASF pools. Table 4 also uses corrected raw-baseline validity checks.
See [run lineage and checks](plinder22-data/minimized-runs/README.md).

## Completed energy and sampling figures, 26 October 2026

Figures 3, 4, and B1 now show CASF core and ref in the approved paired layout.
The artwork uses completed common-protocol energy calculations and exact sampling
curves; review-only preview labels are removed. Results, Methods, captions, and
Supporting Table S9 describe both panels. Table S9 stacks core/ref sections to
keep the energy summaries readable. See [inputs, validation, and reproduction](results-data/core-ref-protocol/README.md).
This supersedes earlier core-only descriptions. The historical compiled previews
are unchanged; review the manuscript in Overleaf.

## Approved tables and Figure 5, 22 October 2026

Tables 1–4 use the refined wide layout (option A), preserving every value and
method order. Table 4 groups reference agreement and pool checks; supplied count
is unranked. Figure 5 now shows PLINDER first, then CASF core/ref, with separate
size and flexibility panels. See [Figure 5 data and reproduction](results-data/size-flexibility/README.md).
The local manuscript preview remains historical; review the compiled update in Overleaf.

## Current PLINDER release, 21 October 2026

All active PLINDER tables and figures now use the cleaned 22-molecule,
1,067-reference release. Figure 1 combines PLINDER diversity/proximity and
coverage/precision; Figure 2 compares CASF core/ref. CASF measurements and
Table 3 are unchanged. See [PLINDER-22 data and reproduction](plinder22-data/README.md).
Earlier figure scripts and provenance are retained as history.

The ACS template now includes a working abstract, Introduction, Methods, Results,
and Discussion, with four main tables, seven main figures, six supporting tables, and six supporting figures.
The separate sampling display is Figure B1. The release checklist (former Table 5) remains removed.
The current combined PLINDER figure is rendered from `plinder22-data/`;
`plinder-threshold-data/` preserves the previous 23-molecule threshold figure.
Title, authors, and end matter still contain upstream examples; Conclusions
remain to be written. This is a working manuscript, not a
submission-ready article.

The former 18-figure review gallery is excluded. The separate Supporting Information
retains the dataset diagram and Qwen scatter, plus four selected analyses refreshed
from core94, ref1044, PLINDER22, the independent torsion runs, and current Qwen
checkpoints. See [appendix-data/README.md](appendix-data/README.md) for selection,
source records, and reproduction. Duplicate supporting tables were removed;
the SI now contains Tables S1–S6 and Figures S1–S6.


## Files to edit

- `introduction.tex`: rationale, relevant literature, and study questions.
- `methods.tex`: benchmark procedures and author comments for missing details.
- `results.tex`: geometric diversity, single-reference recovery, energy dispersion,
  candidate targets, size/flexibility, multiple references, and dataset plans.
- `discussion.tex`: the working Discussion transferred from the Markdown draft.
- `tables/generation-methods.tex`: generator descriptions (Supporting Table S5).
- `tables/diversity.tex`, `tables/count-recovery.tex`, `tables/recovery.tex`,
  `tables/druglike.tex`: Tables 1–4. The unused `tables/release-status.tex` is historical.
- `supporting-information.tex`: current separate SI document with six tables and six figures; see `appendix-data/README.md` for the content map.
- `figures.tex`: figure definitions, captions, and references; artwork is in `figures/`.
- `dataset-figure.tex`: Supporting Figure S1, the experimental-panel diagram, with PDF, SVG, and PNG assets.
- `plinder23-appendix.tex`: selection account and Supporting Tables S2--S4.
- `method-data/`: Methods/appendix drafts, source CSVs, and provenance.
- `results-data/`: current Results/Discussion snapshot, source CSVs, and figure provenance.
- `acs-template.bib`: manuscript references.
- `WRITING_GUIDELINES.md`, `MANUSCRIPT_OUTLINE.md`: agreed writing instructions
  and structure; `METHODS_NOTES.md`, `RESULTS_NOTES.md`: sources and follow-up.

Search for `% AUTHOR` in the LaTeX files for unresolved decisions and analyses.
The notes are included in the Overleaf ZIP but do not appear in the article.

## GitHub and Overleaf workflow

This directory is the manuscript source in
[YerevaNN/casf-benchmark-paper](https://github.com/YerevaNN/casf-benchmark-paper),
on branch `main`. Use the repository-root `main.tex` as the main document in
the linked Overleaf project; retain pdfLaTeX and the Biber bibliography setup.

Future paper edits should be committed and pushed to this repository. Fetch
incoming changes and preserve edits from Overleaf or other authors before
pushing. The older `casf-benchmark/manuscript/` directory is a historical copy.

The current Results artwork, in figure order, is `figure-plinder-combined.*`,
`figure-casf-diversity-recovery.*`, `figure-energy-dispersion.*`,
`figure-3-sampling-budget.*`, and `figure-4-size-flexibility.*` under `figures/`.
PDF, editable SVG, and PNG versions are provided. Main figure numbers are 1–5;
stable filenames need not match display numbers. Current PLINDER and CASF scatter
renderers are listed in `plinder22-data/README.md`. Supporting Figure S1 is the
updated experimental-panel diagram. Earlier figure sources in `results-data/`
and the benchmark `docs/results_figures/` remain historical.

The older `build_figures.py`, `figure-data/`, and `figures/main-figures.pdf`
remain historical assets. They do not reproduce the newly transferred main
figures. Numerical source paths beginning with `docs/`, `data/results/`, or
`src/` refer to the separate `YerevaNN/casf-benchmark` analysis repository.
The 6 October energy addition recalculates MMFF94s energies on the archived
ChEMBL-count ensembles using a common hydrogen preparation without optimizing any coordinates.
It does not regenerate conformers or change the clustering and RMSD results.
The exact protocol, source digests, calculations, and energy renderer are archived
in `energy-data/`. Run `docs/energy_analysis/build_figure.py` in the analysis
checkout after the main figure renderer to update the energy display.

## Preview and verification

`preview/manuscript-preview.pdf` is a historical compiled reading copy of earlier
sections, tables, actual figures, and references. It excludes the
upstream front/end matter. The preview uses Tectonic with biblatex's BibTeX
backend; the Overleaf main document retains Biber, so pagination may differ.
The numerical tables are checked against the corrected exports; size groups
have a reproducible read-only export dated 28 September 2026.

## Template provenance

- [ACS template on Overleaf](https://www.overleaf.com/latex/templates/latex-template-for-american-chemical-society-acs-journal-submissions/swszwgfqsshj)
- [Upstream repository](https://github.com/josephwright/acs-template), retrieved
  25 September 2026 at commit `f9e4afec3eda5d92d7582f459ce0806035d9ebb2`.
- `CC0.txt` and `UPSTREAM_README.md` are unchanged upstream files.
- The section drafts, references, tables, figure slots, and author guidance are
  local additions. The template supports ACS submission, not published layout.

One combined sampling figure (B1) is included after Figure 3, with Hit@0.75
and the RMSD threshold for 80% expected recovery for Qwen, LoQI, FlowR, and
NExT-Mol. Source records, exact subset calculations and the renderer are in
`sampling-data/`; `sampling-figure-options.tex` supplies their captions and method.
Main figure numbering is preserved.

When uploading only the contents of this directory, its `main.tex` also works.
Both entry files load `acs-template.tex`; manuscript, image, and bibliography
paths resolve in either layout.
