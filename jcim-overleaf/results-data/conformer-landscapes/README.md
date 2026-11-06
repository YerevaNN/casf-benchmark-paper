# Conformational landscape visualization

This archive supports the visualization contribution only: main Figure 6 and
its Methods/Results explanation. Existing figure artwork and primary scores retain their content. The molecule-tracking
and Qwen configuration figures are now Figures 7 and 8; Supporting Figures S1–S6
and all tables retain their numbering.
Changes remain uncommitted at the author's request.

## What is presented

The pooled landmark-MDS maps show staurosporine and actinonin for LoQI, FlowR,
and the primary Qwen 1.7B FSQ step-41,054 checkpoint. These examples contrast
concentrated and broader experimental reference collections. All method panels
for a molecule use the same frame and axis limits. The frame was fitted using
all eleven sources from the interactive companion; hiding sources does not
refit it. Up to 70 conformers per source supply the landmarks.

Grey diamonds are cleaned experimental references. Method-specific shading uses
the nearest-reference RMSD, rounded by the page's threshold-preserving export.
Coverage and precision annotations use the existing unrounded manuscript
minima and inclusive 0.75 Å cutoff, independently of the projection.

The only additional geometric summary is the reference cluster count displayed
on these examples: greedy fixed-representative clustering at 0.5 Å, averaged
over 20 random orders. Seeds are 42 + the molecule's alphabetical ordinal in
the full 22-molecule panel. Index 0 is the stored-order control and is excluded
from the displayed mean. Cluster counts describe geometric groupings, not
biological states or energy basins. No reference reweighting, shared-gap
analysis, sampling-budget curves, or new generator ranking is presented here.

## Reproduction

With NumPy and Matplotlib available, run from `jcim-overleaf/`:

```bash
python results-data/conformer-landscapes/build_figure.py
python results-data/conformer-landscapes/verify.py
```

The renderer changes only `figures/figure-conformer-landscapes.{pdf,png,svg}`.
The verifier checks coordinates/point-array alignment, canonical annotation
inputs, and all 42 reference clustering runs. Update output hashes after any
reviewed rendering change.

`landscape-data.json` holds the compact coordinates and nearest-reference
color values for the presented sources. `plot-data.json` records reference
cluster counts and random-order ordinals. `reference_distances.npz` contains
unrounded pairwise RMSD only for the two displayed reference ensembles.
`provenance.json` records input/output hashes, exact model identity, and source
reference-SDF/dataset hashes. Existing `plinder22-data/` remains the numerical
source of all COV-R/P annotations. No conformers were generated or optimized
for this visualization.

The earlier broader analysis and draft integration were retained outside the
paper under the benchmark's `docs/figures/conformational_sampling_story/`.
They are not inputs to the manuscript or the visualization archive.

The figure is rendered at its final 6.5-inch manuscript width, with 7.5-point
metric annotations and generated-pool counts. `landmark-sources.json` identifies
all eleven sources determining the pooled frame. The evaluation appendix links
the interactive companion and records that frames/scales differ across molecules.
