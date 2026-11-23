# CDPKit pharmacophore mixture-policy archive

This directory contains the compact evidence used for the manuscript's final
generator-allocation experiment. The source is the completed, zero-failure
`cdpkit_plinder22_exploration_v1` analysis in the benchmark repository.

The study tests four query families (individual ligand, individual interaction,
shared ligand, and shared protein), seven generation recipes, four attempt
budgets, seven feature requirements, and three positional tolerances. The main
manuscript uses the prespecified stringent six-of-six, 0.5 Å setting. The
allocation decision is based on the individual- and shared-ligand families because
the generators receive no protein structure; the protein-informed families are
secondary validation. The complete primary coverage grid is archived in
`coverage-primary.csv`; it is not reduced to cells that favor a particular
allocation.

Reproduce the archive and manuscript assets from the paper repository root:

```bash
python jcim-overleaf/results-data/cdpkit-pharmacophore-mixture/archive.py
python jcim-overleaf/results-data/cdpkit-pharmacophore-mixture/render.py
```

The rendering step creates main Figure 7, main Table 5, Supporting Tables S9--S10,
and Supporting Figures S11--S12. `provenance.json` records source hashes and
`rendered-files.json` records the compact inputs and rendered outputs.

## Interpretation boundary

The 50:50 allocation is recommended as a balanced default, not as a universal
empirical winner. LoQI leads Qwen at budget 20 in both ligand-only families, while
Qwen gains relative to LoQI as the library expands. The 50:50 allocation exceeds
either pure source in all eight ligand-only family--budget cells, although six
of its eight increments over the stronger pure source have intervals that include
zero. The recommendation combines this crossover with equal source representation
and avoidance of fitting a mixture weight on only 19--21 molecules. Unequal fixed
allocations and the protein-informed families remain in the Supporting Information.
No combined coverage--diversity or coverage--runtime score is constructed, and
generation time is excluded because the generators were not benchmarked under a
common hardware and execution protocol.
