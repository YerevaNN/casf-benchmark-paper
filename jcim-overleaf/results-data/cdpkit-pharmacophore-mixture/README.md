# CDPKit pharmacophore mixture-policy archive

This directory contains the compact evidence used for the manuscript's final
generator-allocation experiment. The source is the completed, zero-failure
`cdpkit_plinder22_exploration_v1` analysis in the benchmark repository.

The study tests four query families (individual ligand, individual interaction,
shared ligand, and shared protein), seven generation recipes, four attempt
budgets, seven feature requirements, and three positional tolerances. The main
manuscript uses the predeclared four-of-four, 0.5 Å setting. The complete primary
coverage grid is archived in `coverage-primary.csv`; it is not reduced to cells
that favor a particular allocation.

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
empirical winner. At budget 20 it is the numerical leader in all four query
families, but its differences from 75:25 are only 0.02--0.19 percentage point and
their molecule-bootstrap intervals include zero. At larger budgets, Qwen or the
75:25 allocation can lead. The recommendation therefore combines the observed
small-budget performance with equal source representation and the avoidance of
tuning a mixture weight on 22 molecules. No combined coverage/diversity score is
constructed.
