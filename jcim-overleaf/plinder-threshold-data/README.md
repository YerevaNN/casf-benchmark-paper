# PLINDER threshold curves (working Figure 5)

The six methods are compared using the existing supplied ensembles for all 23
molecules. A match uses strict RMSD < threshold. Both panels average molecule
fractions equally, without a common PoseBusters filter. No new generation,
optimization, or structural alignment was performed.

- `nearest_distances.npz` and `entries.json`: reference-to-generated and
  generated-to-reference minimum distances for each molecular identity and method.
  These suffice to reproduce any threshold curve. All minima are finite.
- `threshold_curves.csv`: plotted values, 0–1.5 Å at 0.005 Å increments.
- `cache_table_comparison.csv`: per-molecule comparisons at 0.75 Å.
- `provenance.json`: hashes of the original distance matrices and calculation settings.
- `source-preview.py`: original cache extraction/preview script (cluster paths).
- `plot.py`: standalone renderer; run beside the archived CSV and palette.
- `palette.json`: shared manuscript method colors.

Reference coverage is the fraction of reference minima below threshold; precision
is the fraction of generated minima below threshold. Compute each within a molecule,
then average across the 23 molecules. Original matrices contain 2,450 references
per method, with generated counts matching the table records.

The figure is a working descriptive analysis. Cache-derived 0.75 Å summaries differ
from Table 4 by at most 0.0131 percentage points in precision and 0.0046 in coverage.
Table 4 is not overwritten; cache lineage and discrepancies remain to be reconciled
before final evaluation. Common validity filtering is still pending.
