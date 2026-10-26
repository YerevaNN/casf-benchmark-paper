# Figure 6: individual PLINDER molecules across methods

The approved 22-molecule preview is included after Table 4 as Figure 6.
The displayed methods are the previous eight-method selection plus the new
Qwen 0.6B CASF-dedup run. Each line tracks one molecular identity, with the
same color in MAT-R and MAT-P. The three colored molecules are the union of
per-method maxima across the two metrics; the remaining 19 are grey.

`molecule_matching_values.csv` contains all 198 plotted method–molecule rows.
`methods.json` preserves exact identities and plotting order, sorted by mean
MAT-R. The new Qwen label is bold. Input source lineage is preserved in
`source-provenance.json`; `paper-artifacts.json` records the archived inputs and
approved PDF/SVG/PNG copied into `figures/figure-plinder-molecule-tracking.*`.

`plot.py` renders the figure and an interactive companion using these local CSVs.
It writes into this directory; to intentionally refresh the paper, copy its
`plinder22-molecule-tracking.{pdf,svg,png}` outputs to the three figure paths.
Dependencies: NumPy, pandas, matplotlib, Plotly. No alignments or generation are
repeated. The current change copies the already approved preview without rendering
or compiling a manuscript preview.
