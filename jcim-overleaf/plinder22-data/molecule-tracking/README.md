# Figure 6: individual PLINDER molecules across methods

## Primary Qwen checkpoint update, 4 November 2026

All current primary Qwen results use `qwen_1p7b_v2_coordtoken_step41054`, displayed
as **Qwen 1.7B FSQ**. It replaces step47023 in core94, ref1044 and PLINDER22;
the separate Qwen 0.6B CASF-dedup result remains. Tables 1–4, S1–S4/S8/S9,
Figures 1–6/B1 and the optional review figures use the new values. Rankings
exclude the stored ChEMBL3D baseline as before. LoQI now leads PLINDER COV-R
at 0.75 Å; Qwen retains the lowest MAT-R. Training exclusions still require
an executed manifest and must not be inferred from the run name.

Energy and sampling were not published for this checkpoint, so they were measured
from its saved PB-passing conformers under the existing paper protocols. No
coordinate optimization or validity replay was performed. The current energy
totals are 78,058/78,058 core and 734,951/734,952 ref; the one undefined ref
energy is still FlowR's. Ref has 43 empty method–entry energy pools.

See `results-data/qwen-v2/README.md` for exact source identities, reproduction
and validation. This update supersedes older checkpoint statements below;
dated records retain their original lineage. Do not restore old Qwen values
by rerunning a historical exporter.


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
