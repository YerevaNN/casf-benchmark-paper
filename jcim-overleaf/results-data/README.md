# Working Results and Discussion sources

Transferred 5 October 2026 from the benchmark working draft. Only the opening
Results paragraph was substantively revised; scientific values and narrative
order were preserved. Markdown author instructions become LaTeX comments.

- `results-working.md`: source text, table rows, and figure captions.
- `clustering_radius_comparison.csv`, `clustering-radius-provenance.json`: Table 1
  counts at 0.5/1.0 Å, largest-cluster fraction at 1.0 Å, and the read-only
  extraction query and database hash. Generators are sorted by increasing
  0.5 Å count; the stored reference is last.
- `casf_selected.csv`: selected CASF summaries for both candidate targets.
- `druglike_selected.csv`, `druglike_means.csv`: molecule records and summaries.
- `strata.csv`, `paired_reference_coverage.csv`: values underlying Figures 4–5.
- `tables-provenance.json`, `figures-provenance.json`: original input hashes.
- `transfer-provenance.json`: checksums of the transferred Markdown, CSVs, PDFs.
- `source-build-figures.py`: unchanged original renderer, archived for provenance.
  Run the original in the analysis checkout; it assumes that repository layout
  and writes its working Markdown. It is not a paper-repository entrypoint.

Artwork is stored as PDF, SVG, and PNG under `figures/figure-2-*` through
`figure-5-*`. No generation, scoring, or force-field analysis was rerun.
Qwen remains the initial 1.7B FSQ step-47023 checkpoint. Energy results,
benchmark-excluded evaluation, and the public dataset release remain pending.
