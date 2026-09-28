# Manuscript authoring

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

Read METHODS_NOTES.md before revising Methods: it records evidence and unresolved
metadata. methods.tex contains the benchmark Methods draft and is included by
acs-template.tex. Read RESULTS_NOTES.md before revising results.tex: the HTML
report supplies the interpretation, while corrected exports supply the numbers.
introduction.tex is the Introduction draft. Tables are in tables/ and numbered
figure definitions are in figures.tex, with artwork in figures/.
The reproducible renderer is build_figures.py; figure-data/ records the public
release and the denominator differences from dashboard summaries. Methods must describe only the analyses used in Results, tables, or figures.
Energy-window and random-K experiments are outside the current manuscript scope. Front and end matter still contain
upstream examples; Discussion and Conclusions
are not drafted. Replace
sample names and claims only with verified information or explicit placeholders.
For manuscript-only edits, check citations, cross-references, and compilation
when a TeX toolchain is available; Python pipeline tests are not needed.
