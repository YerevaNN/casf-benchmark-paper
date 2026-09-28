# Rendered manuscript preview

`manuscript-preview.pdf` is a 11-page reading copy of the drafted Introduction,
Methods, and Results, with three tables, four figures, and nine
references. It imports the section files directly and excludes the upstream
template's example front/end matter. Page renders may be generated locally for review.

Compiled locally with Tectonic 0.17.0. The preview uses biblatex's BibTeX
backend and a small line-breaking allowance; the Overleaf main document retains
Biber. Pagination may differ. The rendered PDF was checked for all table/figure
numbers, bibliography entries, and unresolved-reference markers. The BibTeX
fallback still emits a generic bibliography rerun warning, although the
rendered citations and manuscript cross-references resolve.

Run `tectonic manuscript-preview.tex` from this directory to refresh when
Tectonic is available. Generated previews are excluded from the Overleaf ZIP.
