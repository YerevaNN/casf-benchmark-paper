# Selected Supporting Information

Table S8 now records completed v2 Qwen training exposure; detailed training and
evaluation appendices follow S7. Both PDF entry points include all eight tables and
six figures. Training evidence is in `../method-data/training-audit/`.

The former 18-figure review gallery is no longer included in the manuscript.
The current SI contains seven tables and six figures, all linked to the main text.
It is appended to the main review PDF and also available as a standalone document;
both entry points include `supporting-content.tex`.
Existing measurements are reorganized and summarized here; no conformers were
regenerated, optimized, filtered, or realigned.

## Retained figures

| Current figure | Purpose | Previous material |
| --- | --- | --- |
| S1 | Evaluation panels and identity selection | Existing dataset diagram |
| S2 | Per-molecule Qwen representation and size effects | Existing approved scatter |
| S3 | Recoveries gained and lost versus ChEMBL3D | C02, expanded to core94 and ref1044 |
| S4 | Joint size/flexibility composition and sparse groups | C07, expanded to both CASF panels |
| S5 | Per-molecule coverage and precision hidden by the averages | C11, expanded to all eleven Table 4 methods |
| S6 | Per-molecule validity and the four-Qwen comparison | C14, expanded to Table 4 methods and all four v2 configurations |

S3 uses the fixed 1,000-candidate pools after PB and the unchanged stored baseline.
Missing hits count as misses; all 94/1,044 identities remain in the denominator.
S4 cross-tabulates the same reference descriptors used in main Figure 5; its
margins are checked against that figure. S5 uses inclusive RMSD <= 0.75 Å and the
full supplied PLINDER pools. S6 measures validity separately. Heatmap values are
rounded to whole percentages, while all calculations retain full precision.

The two PLINDER torsion rows are `torsion_independent_raw_druglike` and
`torsion_independent_minimized_druglike`. Primary Qwen is
`qwen_1p7b_v2_coordtoken_step41054`. The older 0.6B CASF-dedup step64000 model
is distinguished from the v2 binned step66057 configuration. The four-v2 source
archive records the precise labels and generation/training metadata.

## Tables and removals

Table S1 now gives retained counts for both CASF panels and both candidate
targets; recovery percentages already appear in main Table 3. Tables S2–S4
retain PLINDER selection criteria, statistics, and identities. Table S5 retains
generator provenance and S6 retains the common-protocol energy summaries.
Table S7 adds the verified core/ref selection counts; its source logs are in
`method-data/selection-audit/`. Existing table and figure numbering is unchanged.

Former S1 (core minimum/median RMSD), S3 (size-stratified recovery), and S4
(PLINDER matching means) are excluded: the first adds an unused summary, and
the latter two repeat the main size figure and Table 4. Former S2 is replaced
by the retained-count table. Original source files remain historical records.

C01/C17 bootstrap alternatives, C03/C04 recovery/complementarity matrices,
C05/C06 alternative distance comparisons, C08 diversity alternatives, C09
aggregate validity, C10 separately sampled classical endpoint plots, C12/C15
aggregate PLINDER plots, C13/C18 selected molecule comparisons, and C16
reference counts are excluded. Their information is either already in the main
text/current SI, less directly related to its claims, or redundant with the
retained full molecule heatmaps. The old assets and renderer are archived;
they are not inputs to the selected figures or compiled documents.

## Reproduction

With the benchmark environment activated, run from `jcim-overleaf/`:

```bash
python appendix-data/build.py
python appendix-data/verify.py
python plinder22-data/verify.py
```

`build.py` reads only the current manuscript CSV archives. It validates their
aggregates against Tables 3/4, Figure 5 strata, and Figure 7, then creates the
selected PDF/PNG/SVG artwork, supporting counts table, and exact plotted CSVs.
`provenance.json` records input/output hashes. Main-text values and figures are
not rebuilt. Manuscript compilation remains an Overleaf review step.
