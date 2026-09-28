# Manuscript figures

Four figures are included in the paper, with captions and placement controlled
by `../figures.tex` and `../results.tex`:

1. `recovery-thresholds-budget`: RMSD thresholds and candidate-budget endpoints.
2. `diversity-recovery`: cluster radii of 0.5, 1.0, and 2.0 Å and diversity versus recovery.
   The archived data retain the 3.0 Å results, which are omitted from the plot.
3. `size-flexibility`: core/ref molecular size and flexibility groups.
4. `multireference-coverage`: drug-panel recall/precision with molecule bootstrap
   intervals and contrasting paired examples.

Each has a vector PDF, editable SVG, and 450-dpi PNG at 7 inches wide.
The paper includes the PDFs. Sources and numerical exports are in
`../figure-data/`; `../build_figures.py` reproduces the outputs.

`energy-window-historical` is a separate diagnostic of the older extended
sidecar in the public release. Its legacy Qwen identity and historical status
are explicit; it is NOT included in the manuscript or the main figure gallery.
Current-checkpoint energy windows and fixed-valid-K curves still need updating.

Formatting follows the size, font-legibility, and accessible-color guidance in
[JCIM's graphics instructions](https://researcher-resources.acs.org/publish/author_guidelines/pdf?coden=jcisd8).
