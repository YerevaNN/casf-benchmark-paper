# JCIM manuscript

The ACS template now includes drafted Introduction, Methods, and Results
sections, three tables, and four figures. Title, authors,
abstract, and end matter still contain upstream examples. Discussion and
Conclusions remain to be written. This is a working manuscript, not a
submission-ready article.

## Files to edit

- `introduction.tex`: rationale, relevant literature, and study questions.
- `methods.tex`: benchmark procedures and author comments for missing details.
- `results.tex`: recovery-first story using core94 and Qwen 1.7B FSQ step47023.
- `tables/recovery.tex`, `tables/size.tex`, `tables/druglike.tex`: included tables.
- `figures.tex`: figure definitions, captions, and references; artwork is in `figures/`.
- `acs-template.bib`: nine verified references.
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

The four figure PDFs are included, with a combined gallery in
`figures/main-figures.pdf`. Each also has an editable SVG and a 450-dpi PNG.
Run `python jcim-overleaf/build_figures.py` from the repository root to render
from the archived CSVs; this does not require the benchmark databases.
The script's `--refresh-data` mode was written for the benchmark checkout and
requires access to its source databases; do not refresh the evidence silently.

Numerical source paths beginning with `docs/`, `data/results/` or `src/` in
notes refer to the separate `YerevaNN/casf-benchmark` analysis repository.
Figure-source CSVs and the public release record are archived here under
`figure-data/`. Figure 1 uses candidate-tier endpoints rather than historical
random-K curves; the old energy analysis remains a separate diagnostic.

## Preview and verification

`preview/manuscript-preview.pdf` is a compiled reading copy of the authored
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
