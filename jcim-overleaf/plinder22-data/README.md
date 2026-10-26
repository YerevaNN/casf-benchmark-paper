# Cleaned PLINDER-22 paper results

The 20 October 2026 evaluation replaces PLINDER-23 in all active manuscript
PLINDER displays. The original eight methods retain their generated pools, and two minimized
classical variants are added. Qwen retains its 1.7B FSQ step47023 identity;
no retrained-model result is substituted in that ten-method archive.
Table 4 and S4 additionally include the Qwen 0.6B CASF-dedup run archived
in `qwen-dedup/`; Figure 6 uses the approved molecule-tracking subset
in `molecule-tracking/`. The new run is not added to the other displays.

The cleaned panel excludes bortezomib and three wrong-stereo reference records
(two chloramphenicol and one staurosporine). It contains 1,067 representatives
across 22 molecules, 552 PDB entries, and 1,045 systems. Representatives were
supplied after 0.05 Å deduplication and anchor-triplet alignment; duplicate-cluster
sizes are not evaluation weights. Full upstream cleaning code/release lineage
remains outstanding. `reference-source-notes.md` preserves the delivery audit.

See [finalized minimized-run provenance](minimized-runs/README.md) for the paired
minimization protocol and corrected raw-baseline validity records.

Figure 1 uses A across the top (all ten methods) and B/C side by side below
(RDKit raw plus the six learned methods). All ten methods remain in the archived
curves and tables; the figure restriction is a display choice.

## Data and definitions

- `per-molecule.csv` and `summary.csv`: selected rows from the refreshed
  benchmark `data/results/plinder_panel/` exports. Cluster counts and occupancy
  use full supplied pools. PB validity is separate; no PB mask is applied to
  clustering or COV/MAT. Table 4 pools PB counts across molecules.
- `nearest-distances.npz`, `entries.json`, `threshold-curves.csv`: refreshed
  reference/generated minima and equal-molecule-weight curves. A match is
  RMSD **<=** threshold. Every reference and generated minimum is finite.
  Tables 2, 4, S4 and all current PLINDER curves use these same measurements.
- `reference-inventory.csv` and `reference-records.csv`: counts and complete SDF
  record properties, including PDB/system IDs, original indices and cluster
  membership. SDF hashes are in `provenance.json`; structures are not duplicated.
- `methods.json`: exact method identities and display names.
- `provenance.json`: source paths, SHA-256 hashes, evaluation conventions and
  checks against all ten refreshed caches. The older `druglike_cleaned` tables
  are not mixed in: the panel export supplies full-pool clustering and occupancy.
- `rendered-files.json`: hashes of current display inputs, renderers and outputs.

## Reproduction

From the paper repository, with NumPy/pandas/Matplotlib available:

```bash
python jcim-overleaf/results-data/three-dataset-tables/build_tables.py
python jcim-overleaf/results-data/build_druglike_table.py
python jcim-overleaf/plinder22-data/build_inventory.py
python jcim-overleaf/plinder22-data/plot.py
python jcim-overleaf/plinder22-data/plot_datasets.py
python jcim-overleaf/results-data/three-dataset-tables/plot_casf.py
python jcim-overleaf/build_candidate_figures.py --plinder-only
python jcim-overleaf/plinder22-data/verify.py
```

`export.py /path/to/casf-benchmark --minimized-benchmark /path/to/minimized-checkout
--baseline-runs /path/to/druglike_baselines` refreshes the archived inputs from the local
benchmark release and requires its environment/RDKit. It validates cache shapes,
identities, finite minima, COV/MAT at each headline cutoff, the published curves,
and cluster summaries before replacing the PLINDER table data. It does not run
new generation, alignment, clustering or energy calculations. CASF CSV cells
are preserved.

The combined PLINDER figure is main Figure 1; CASF core/ref diversity–proximity
is Figure 2. Energy and sampling figures compare CASF core/ref; Figure 5 includes all three datasets
and all ten PLINDER methods.
Table 3 remains CASF-only. PLINDER review candidates C11–C18 retain the six learned
methods; C18 keeps the three previously chosen examples that remain in the panel.
The earlier PLINDER-23 figures and dated audit files are historical and are not
included as current results. Tables S6/S7 keep their original `plinder23` filenames
and labels for stable includes; the current identity CSV is
`../method-data/plinder22-molecules.csv`.

Validation covers numeric consistency, table rendering, source includes, labels,
citations, artwork existence and visual inspection. No full manuscript was compiled
locally, following the paper workflow; layout in the compiled paper is reviewed in Overleaf.
