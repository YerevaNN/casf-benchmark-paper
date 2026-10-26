# JCIM manuscript

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
and Discussion, with four main tables, five main figures, nine supporting tables, and one supporting figure.
The separate sampling display is Figure B1. The release checklist (former Table 5) remains removed.
The current combined PLINDER figure is rendered from `plinder22-data/`;
`plinder-threshold-data/` preserves the previous 23-molecule threshold figure.
Title, authors, and end matter still contain upstream examples; Conclusions
remain to be written. This is a working manuscript, not a
submission-ready article.

The draft also includes **18 additional figure candidates (C01–C18)** in a
review appendix. See [FIGURE_CANDIDATES.md](FIGURE_CANDIDATES.md) for the index
and [the gallery](figures/candidates/additional-figures.pdf) for all options.
They visualize the HTML report's questions using corrected archived measurements;
the candidates remain separate from the current main figures. Each candidate has a full caption in
`figure-candidates.tex`. Remove individual blocks after selection, or comment out
its input in `acs-template.tex` to hide the whole review appendix.
Recreate this set with `python jcim-overleaf/build_candidate_figures.py`.

## Files to edit

- `introduction.tex`: rationale, relevant literature, and study questions.
- `methods.tex`: benchmark procedures and author comments for missing details.
- `results.tex`: geometric diversity, single-reference recovery, energy dispersion,
  candidate targets, size/flexibility, multiple references, and dataset plans.
- `discussion.tex`: the working Discussion transferred from the Markdown draft.
- `tables/generation-methods.tex`: generator descriptions (Supporting Table S8).
- `tables/diversity.tex`, `tables/count-recovery.tex`, `tables/recovery.tex`,
  `tables/druglike.tex`: Tables 1–4. The unused `tables/release-status.tex` is historical.
- `supporting-information.tex`: separate SI document with Tables S1--S9 and Figure S1 (core RMSD,
  ref recovery, size strata, multi-reference matching distances, and the PLINDER
  selection criteria, statistics, molecule list, generator descriptions, and energy summaries); switch the Overleaf main document to this file to compile it.
- `figures.tex`: figure definitions, captions, and references; artwork is in `figures/`.
- `dataset-figure.tex`: Supporting Figure S1, the experimental-panel diagram, with PDF, SVG, and PNG assets.
- `plinder23-appendix.tex`: selection account and Supporting Tables S5--S7.
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
