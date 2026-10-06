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


The author's latest scope decision is core94 and a single recovery cutoff of
0.75 Å. Keep ref in SI and remove threshold-sensitivity narratives from the
main paper. Emphasize selection among the broader method panel under the
study's identity, validity, and candidate-budget constraints. Table 1 reports 0.5 and 1.0 Å
clustering radii, sorted by the 0.5 Å mean; other diversity displays use 1.0 Å.
These remain distinct from the recovery cutoff. See the updated
MANUSCRIPT_OUTLINE.md for current figure and table placement.

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
Energy-window and random-K experiments are outside the current manuscript scope. Front and end matter still contain
upstream examples; Discussion is a working draft and Conclusions remain
to be written. Replace
sample names and claims only with verified information or explicit placeholders.
For manuscript-only edits, check citations and cross-references in the source.
Do not compile local previews: the user reviews changes in Overleaf.
Python pipeline tests are not needed.
