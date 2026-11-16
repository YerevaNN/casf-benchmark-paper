# Figure 7: individual PLINDER molecules across methods

The figure follows 22 molecules across eight generators: Qwen 1.7B FSQ,
five external learned methods, and two raw classical baselines. The primary
Qwen model is `qwen_1p7b_v2_coordtoken_step41054`.

`molecule_matching_values.csv` contains 176 method–molecule rows. `methods.json`
records exact identities and order by mean MAT-R. Lines connect the same molecule;
the union of per-method maxima across MAT-R/P is highlighted. All values come
from the supplied ensembles, without a common PoseBusters filter.

Run `plot.py` with the benchmark environment, then copy its PDF/PNG/SVG outputs
to `figures/figure-plinder-molecule-tracking.*`. The interactive HTML is a local
companion. No generation, minimization, alignment or scoring is repeated.
`paper-artifacts.json` records current inputs and figure hashes;
`source-provenance.json` preserves the original upstream extraction lineage.
