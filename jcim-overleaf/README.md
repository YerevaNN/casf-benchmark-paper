# JCIM manuscript

The ACS template now includes a working abstract, Introduction, Methods, Results,
and Discussion, with five main tables, five main figures, nine supporting tables, and one supporting figure.
Table 5 records pending dataset-release items and is not a measured result.
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
  `tables/druglike.tex`, `tables/release-status.tex`: Tables 1–5.
- `supporting-information.tex`: separate SI document with Tables S1--S9 and Figure S1 (core RMSD,
  ref recovery, size strata, multi-reference matching distances, and the Plinder-23
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
on branch `main`. Use `jcim-overleaf/acs-template.tex` as the main document in
the linked Overleaf project; retain pdfLaTeX and the Biber bibliography setup.

Future paper edits should be committed and pushed to this repository. Fetch
incoming changes and preserve edits from Overleaf or other authors before
pushing. The older `casf-benchmark/manuscript/` directory is a historical copy.

The current Results artwork is `figures/figure-2-diversity-recovery.*`,
`figure-energy-dispersion.*`, `figure-3-sampling-budget.*`, `figure-4-size-flexibility.*`, and
`figure-5-multiple-references.*`, each in PDF, editable SVG, and 350-dpi PNG.
Supporting Figure S1 is the experimental-panel diagram. Main figure numbers are
1–5; the artwork filenames retain their original stable identifiers. The current figures
come from the benchmark's `docs/results_figures/build_figures.py`. The initial
renderer is archived as `results-data/source-build-figures.py`; the 6 October
revision, with Figure 1 comparing 0.5 Å cluster count and Best RMSD, is archived
as `results-data/source-build-figures-2026-10-06.py`.
That script runs in the original analysis checkout and also updates its
working Markdown; it is not a standalone renderer for this paper checkout.

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
