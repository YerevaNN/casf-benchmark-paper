# CDPKit pharmacophore mixture-policy archive

This directory contains the compact evidence used for the manuscript's final
generator-allocation experiment. The source is the completed, zero-failure
`cdpkit_plinder22_exploration_v1` analysis in the benchmark repository.

The study tests four query families (individual ligand, individual interaction,
shared ligand, and shared protein), seven generation recipes, four cached-conformer
budgets, seven feature requirements, and three positional tolerances. Here, a budget
counts accepted cached conformers rather than generation attempts. The main
manuscript uses the prespecified stringent six-of-six, 0.5 Å setting. The
allocation decision is based on the individual- and shared-ligand families because
they represent the intended reusable, target-independent library; the protein-informed
families test binding-site relevance as secondary validation. The complete primary coverage grid is archived in
`coverage-primary.csv`; it is not reduced to cells that favor a particular
allocation.

Qwen's completed cache is separate from this sampling budget. Its original cache
contained 21,947 accepted conformers, and a completion pass made 59 additional
decoding attempts to retain 53 additions. The exact record is archived in
`../pharmacophore-query/pilot-qwen-topup-manifest.json`; corresponding LoQI timing
was not measured under a common execution protocol.

Reproduce the archive and manuscript assets from the paper repository root:

```bash
python jcim-overleaf/results-data/cdpkit-pharmacophore-mixture/archive.py
python jcim-overleaf/results-data/cdpkit-pharmacophore-mixture/render.py
```

The rendering step creates main Figure 7, main Table 5, Supporting Tables S9--S10,
and Supporting Figures S11--S12. `provenance.json` records source hashes and
`rendered-files.json` records the compact inputs and rendered outputs.

## Interpretation boundary

The 50:50 allocation is presented as a reasonable fixed production choice, not as
a uniquely established optimum. LoQI leads Qwen at budget 20 in both ligand-only families, while
Qwen gains relative to LoQI as the library expands. The 50:50 allocation exceeds
either pure source in all eight ligand-only family--budget cells, although six
of its eight increments over the stronger pure source have intervals that include
zero. The recommendation combines this crossover with equal source representation
and a simple equal-allocation rule. The data do not reliably distinguish 50:50 from
75:25 at the intended larger budgets. Diversity is a
secondary selection axis: at budget 1,000, 50:50 has 220.5 mean clusters versus
100.7 for LoQI, 256.7 for 75:25, and 279.6 for Qwen. Thus 50:50 preserves substantial
breadth relative to LoQI but does not maximize diversity. The Qwen-heavy 75:25 allocation
is shown in every main-figure panel; 25:75, the full seven-recipe endpoint table,
and the complete grid remain in the Supporting Information. No combined
coverage--diversity or coverage--runtime score is constructed, and
generation time is excluded because the generators were not benchmarked under a
common hardware and execution protocol.
