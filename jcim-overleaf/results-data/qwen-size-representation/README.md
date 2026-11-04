# Qwen size and representation comparison

The approved aggregate comparison is Figure 7, after the PLINDER analysis.
The simplified molecule-level scatter is Supporting Figure S2 in
`supporting-information.tex`. Earlier tables and Figures 1–6 are unchanged.
PDF and PNG artwork was copied byte-for-byte from the approved previews; SVG
trailing whitespace was normalized without changing the artwork. No conformers
were generated, optimized, filtered or realigned for this manuscript addition.

## Exact models

| Key | Run |
|---|---|
| binned06 | `qwen_0p6b_v2_binned_pm12_step66057` |
| binned17 | `qwen_1p7b_v2_binned_pm12_step66057` |
| fsq06 | `qwen_0p6b_v2_fsq_step41054`, corrected export |
| fsq17 | `qwen_1p7b_v2_coordtoken_step41054` |

The 0.6B binned model is **not** the CASF-dedup PM12 step64000 model in Table 4
and Figure 6. FSQ 1.7B is the primary model used elsewhere in the paper.

`sources.json` preserves upstream paths and hashes from the preview analysis;
its original scope string records that earlier preview stage.
`provenance.json` records the approved manuscript assets and archived inputs.
Training configuration snapshots and source hashes are in
`training-configurations.json`; generation metadata is archived separately.
Within each representation, the sizes share the data manifest and epoch plan.
Across representations, training exports and token exposure differ. These are
comparisons of trained configurations, not a pure causal ablation of encoding.
Executed exclusion manifests remain a separate author follow-up.

## Evaluation and observations

CASF uses 1,044 ref identities at the fixed 1,000-candidate target after PB
filtering. All four configurations have measurable distances for every identity;
retained counts differ. PLINDER uses the same 22 molecules and 1,067 references
for each configuration. Its matching statistics use full supplied pools, with
PB measured separately. Matching averages weight each identity equally; the
displayed PB rate pools passes and supplied counts across molecules.

On ref, FSQ improves matching for 695/1,044 (66.6%) and 687/1,044 (65.8%)
identities at 0.6B and 1.7B. Scaling improves 527/1,044 (50.5%) for binned and
563/1,044 (53.9%) for FSQ. The same underlying values appear in both main-figure
rows; the rows change which configurations are connected, not the observations.

On PLINDER, the pooled FSQ PB pass rate increases from 82.7079% to 92.4773%.
Chloramphenicol and Staurosporine account for 86.0874% of the **net sum of
per-molecule percentage-point changes**, equivalently the net change in their
unweighted mean. This contribution is not a decomposition of the pooled rate.
The underlying values are in `validity-by-molecule.csv` and
`plinder-per-molecule.csv`.

Actinonin contributes to the mean PLINDER representation difference. Including
all molecules, binned-minus-FSQ mean MAT-R is 0.04314 Å at 0.6B and 0.03207 Å
at 1.7B; omitting Actinonin in a sensitivity calculation yields 0.01731 Å and
0.00506 Å. The manuscript figures retain all molecules. No uncertainty intervals
or claims of training-seed reproducibility are added by this comparison.

## Reproduce the artwork

Using the analysis repository's activated environment (pandas, NumPy,
Matplotlib), run:

```bash
python results-data/qwen-size-representation/build_figures.py
```

Run from `jcim-overleaf/`. The wrapper renders from the archived CSVs and copies
PDF, SVG and PNG assets into `figures/`. The CASF archive retains only fixed-target
ref rows and the columns used by these figures. `summary.csv` preserves the
original summary export, but only ref/fixed and PLINDER/supplied are plotted.

Validation: exact approved-asset checksums; matched identities and reference
counts; finite distances; recomputed means, pooled PB rates and paired changes;
source citation/cross-reference and figure-path checks. No local manuscript
compilation was performed; review the main manuscript and separate Supporting
Information in Overleaf.
