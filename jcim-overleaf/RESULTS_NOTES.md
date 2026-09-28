# Results: sources and author follow-up

Revised 28 September 2026. The manuscript now follows the agreed recovery-first
story: core94 anchors the comparison; Qwen 1.7B FSQ pretraining step47023 is
the sole main-text Qwen; ref supports specific findings; the 23-molecule
multi-reference drug panel is a separate experiment. Additional Qwen variants
and SFT comparisons are reserved for possible Supporting Information.

## Interpretation and numerical sources

The biologist's `docs/jcim_publication_report_2026_09_14.html` supplies the
interpretive perspective, especially its distinction between recovery,
sampling efficiency, diversity, and multi-reference precision. Its older
numbers do not supersede the corrected exports.

The source database fingerprint was checked against the 24 September snapshot
before the size export. No generation, energy evaluation, or CASF random-K
analysis was rerun in this manuscript revision.

| Content | Source relative to repository root |
| --- | --- |
| Recovery, threshold crossings, retained counts, minimum/median RMSD, diversity | `docs/publication_tables_2026_09_24/all_casf_summary.csv` |
| Paired intervals and cohort-overlap sensitivity | `docs/publication_tables_2026_09_24/paired_cluster_bootstrap.csv` |
| Smaller sampling targets | `docs/publication_tables_2026_09_24/tiers.csv` |
| Size strata and Table 2 | `docs/publication_tables_2026_09_28/size_strata.csv`, reproduced by `docs/publication_size_analysis_2026_09_28.py` |
| Flexibility | `docs/publication_tables_2026_09_24/flexibility.csv` |
| Table 3 coverage/MAT | Molecule means of `docs/publication_tables_2026_09_24/druglike_per_molecule.csv` |
| Table 3 PB | Pooled pass fraction in `docs/publication_tables_2026_09_24/druglike_summary.csv` |
| Drug uncertainty and examples | `docs/publication_tables_2026_09_24/druglike_paired.csv` and `druglike_cases.csv` |

Table 1 contains the full selected generator panel plus ChEMBL3D-PB, with
core retained counts and RMSD. Table 2 uses four contrasting methods for
readability; the export contains the wider panel. Table 3 includes all five
external methods and the selected Qwen. All three tables are separate LaTeX
inputs under `tables/` and are included by `results.tex`.

A material threshold result is retained: on core, LoQI exceeds Qwen at 0.5 Å
(83.0% versus 78.7%) despite Qwen's higher recovery at 0.75 Å. Therefore the
text does not claim threshold-invariant superiority. No correlation from the
old 19-generator panel is presented as a result for the reduced main panel.

## Figures recreated from the public dashboard

On 28 September, GitHub's latest public dashboard release was verified as
`dashboard-data-stereo-identity-v2` (18 September). Both database digests match
the local files exactly. Figure sources and explicit denominator/checkpoint
choices are archived in `figure-data/README.md` and `provenance.json`.

Four real figures now replace the placeholders and are inserted near their
discussion using commands defined in `figures.tex`:

1. `recovery-thresholds-budget.pdf`: core RMSD thresholds and ChEMBL-count versus
   fixed candidate endpoints. These are not random valid-K curves.
2. `diversity-recovery.pdf`: clustering-radius curves and cluster count versus
   recovery. It does not contain an unverified current energy panel.
3. `size-flexibility.pdf`: core and supporting ref strata, with group sizes.
4. `multireference-coverage.pdf`: raw-pool recall/precision with 95% marginal
   molecule-bootstrap intervals and paired molecule examples.

PDFs, SVGs, and 450-dpi PNGs are in `figures/`. The public release's unchanged
14 September energy sidecar was recreated as `energy-window-historical.pdf`,
clearly labeled historical and excluded from the paper. No current energy
robustness result is inferred from it. Additional Qwen labels in the release
were not substituted for the exact step47023 checkpoint.

The figure script also archives the dashboard's displayed summaries. MCF and
NExT-Mol core Hit@0.75 differs from the displayed conditional percentage
(81/93 = 87.1%) because the paper counts the missing entry as a miss
(81/94 = 86.2%). All manuscript figures retain the paper's defined denominator.

## Remaining scientific decisions

- Energy windows, dynamic tiers, normalized entropy, and drug random-K
  analysis are outside the current story; their Methods and energy-related
  Introduction/Results prose were removed. Historical artifacts remain
  archived but are not manuscript experiments.
- Confidence intervals were calculated from the underlying records: paired
  statistics in the corrected 24 September analysis, and marginal Figure 4
  intervals during figure preparation. See METHODS_NOTES.md for provenance.
- Harmonize drug PB filtering and RMSD-failure handling; update prose, tables,
  and figures together. Current values are descriptive supplied-pool results.
- Document representative-checkpoint selection, training/inference metadata,
  training-set overlap, and dataset provenance (see METHODS_NOTES.md).
- Keep sparse size strata descriptive. Size is not an independent causal
  explanation after stratifying on it alone.
- Additional Qwen/SFT tables, dynamic tiers, and full threshold combinations
  can go in SI; no nonexistent SI table numbers are cited.

The Introduction cites prior experimental-conformation benchmarking rather
than claiming priority. The added McNutt et al. reference was verified against
https://pubs.acs.org/doi/10.1021/acs.jcim.3c01245 on 28 September 2026.
