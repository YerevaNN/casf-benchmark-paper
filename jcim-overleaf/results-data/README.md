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
- `source-build-figures-2026-10-06.py`: revised renderer with Figure 1 showing
  0.5 Å cluster count versus Best RMSD (mean of per-entry minima), and the
  revised Figure 3 using 0.5 Å clusters and Hit@0.75 in Table 3 row order.
- `diversity-rmsd-points.csv`: exact coordinates and measurable-entry counts for
  Figure 1.

Artwork is stored as PDF, SVG, and PNG under `figures/figure-2-*` through
`figure-5-*`. Clustering and RMSD measurements are unchanged. The new energy calculation
and Figure 2 are archived separately in `../energy-data/`.
Qwen remains the initial 1.7B FSQ step-47023 checkpoint. The benchmark-excluded evaluation and public dataset release remain pending.

Table 3 and Figure 3 use `sampling_budget_comparison.csv`, exported by
`build_sampling_budget.py` from the pinned database. `sampling_budget_per_entry.csv`
archives the underlying entry records and `sampling-budget-provenance.json`
records the input and output hashes. Run the script in the analysis checkout.
The 1,000-candidate columns precede ChEMBL-count columns. The stored ChEMBL3D
ensemble is shown only in the latter; no larger reference pool is invented.
