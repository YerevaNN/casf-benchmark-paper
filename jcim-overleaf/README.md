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

## Manual Overleaf upload

1. Upload `jcim-overleaf.zip` as a new Overleaf project.
2. Select `acs-template.tex` as the main document and use pdfLaTeX.
3. Recompile; Overleaf should run Biber for the `biblatex` bibliography.

The four figure PDFs are already included. `figures/main-figures.pdf` provides
a compact gallery of the main figures. Each also has an editable SVG and
a 450-dpi PNG under `figures/`. `build_figures.py` recreates them from the
archived public-dashboard values in `figure-data/`, which records release
identity, database hashes, checkpoint selection, and denominators.

Figure 1 uses corrected candidate-tier endpoints rather than historical
random-K curves. Figure 2 uses current diversity results; the older energy
analysis is saved separately and clearly labeled historical. The drug figure
still describes supplied pools pending common validity/failure handling.

The ZIP is a snapshot. To rebuild it from the repository root:

```python
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
root = Path('manuscript')
with ZipFile(root / 'jcim-overleaf.zip', 'w', ZIP_DEFLATED) as archive:
    for path in sorted(root.rglob('*')):
        if not path.is_file() or 'preview' in path.relative_to(root).parts:
            continue
        if path.suffix in {'.tex', '.bib', '.md', '.txt', '.py', '.csv', '.json'} or 'figures' in path.relative_to(root).parts:
            archive.write(path, path.relative_to(root))
```

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
