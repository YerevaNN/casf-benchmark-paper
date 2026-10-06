# Results

> **Data status:** Tables use the corrected findings archived on 24 September 2026. Qwen values belong to the initial 1.7B FSQ checkpoint at step 47,023; they do not describe the pending evaluation after benchmark exclusions. Energy results use the separate common-hydrogen rescoring; release counts remain pending.

A useful conformational resource should represent a range of molecular shapes while retaining coverage of geometries observed in protein complexes. We begin by examining the variation present in the generated ensembles, before asking whether that variation includes the experimental bound conformations. The first comparison uses the 94-entry CASF core panel at the ChEMBL-count target, placing each generator alongside the stored ChEMBL3D ensemble for the same molecules. Table 1 summarizes geometric diversity through cluster counts at 0.5 and 1.0 Å, together with the retained ensemble sizes and the fraction of conformers in the largest cluster.

<!-- results-table-2:start -->
**Table 1. Geometric diversity at the ChEMBL-count target on the 94-entry core panel.**

| Method | Mean retained conformers ↑ | Mean clusters at 0.5 Å ↑ | Mean clusters at 1.0 Å ↑ | Mean largest-cluster fraction at 1.0 Å (%) ↓ |
| --- | --- | --- | --- | --- |
| RDKit minimized | **81.2** | 21.7 | 11.0 | 49.8 |
| Torsion minimized | **81.2** | 21.9 | 11.7 | 56.4 |
| LoQI | **81.2** | 25.0 | 9.2 | 53.7 |
| MCF drugs-L | 66.1 | 31.2 | 14.2 | 47.3 |
| NExT-Mol DMT-L | 66.4 | 31.6 | 14.2 | 46.4 |
| RDKit raw | **81.2** | 34.3 | 14.7 | 46.1 |
| Torsional Diffusion | 69.3 | 40.9 | 17.5 | 44.7 |
| FlowR | 72.7 | 54.1 | **29.1** | **42.4** |
| Torsion raw | 71.1 | <u>54.2</u> | 21.8 | 44.4 |
| Qwen 1.7B FSQ | <u>78.3</u> | **55.3** | <u>22.7</u> | <u>44.1</u> |
| ┄┄┄ | ┄┄┄ | ┄┄┄ | ┄┄┄ | ┄┄┄ |
| ChEMBL3D-PB | **81.2** | 21.7 | 7.7 | 47.2 |

Generators are ordered by increasing unrounded mean cluster count at 0.5 Å; the stored ChEMBL3D-PB ensemble is shown separately. Bold and underlining mark the best and second-best distinct values in each column, including ties and the stored reference. More clusters indicate broader geometric diversity; a smaller largest-cluster fraction indicates less concentration in one cluster. Retained count describes yield, not diversity. Retained counts are averaged over all 94 entries, including empty outputs. Cluster statistics use identical defined-entry sets within each method: 92 for NExT-Mol, 93 for MCF, and 94 for the others. Largest-cluster fractions are calculated per molecule before averaging. These are geometric clusters, not energy basins. Source: [clustering and occupancy comparison](results_tables/clustering_radius_comparison.csv).
<!-- results-table-2:end -->

Qwen, random torsion sampling, and FlowR produce the most geometrically diverse ensembles by cluster count at both radii. Their counts are close at 0.5 Å, with Qwen ranked first, whereas FlowR produces the most clusters at 1.0 Å. LoQI, the minimized baselines, and the stored ChEMBL3D ensembles consistently form fewer clusters at both radii, indicating narrower geometric variation under these sampling conditions. The broad separation between these groups therefore persists across the two resolutions, despite changes in their exact ordering.

Cluster occupancy provides a complementary description of how the samples are distributed. At 1.0 Å, FlowR and Qwen have the smallest mean fractions of conformers in the largest cluster, at 42.4% and 44.1%, respectively; random torsion sampling is close at 44.4%. LoQI and minimized torsion sampling are more concentrated, with 53.7% and 56.4% in the largest cluster. ChEMBL3D has fewer clusters but a less concentrated largest cluster than these two methods, illustrating that the number of represented shapes and their sampling balance describe different aspects of diversity. These observations concern the geometry of the retained pools, whose sizes differ between methods; they do not establish whether the additional shapes are energetically plausible or close to experimentally observed conformations.

[Update the method-specific observations and their magnitudes from the final ChEMBL-count evaluation; do not transfer rankings from the 1,000-candidate comparison to this target.]

Having characterized the ensembles themselves, we assess recovery of the experimental bound geometries at the ChEMBL-count target. Table 2 reports Hit@0.5 and Hit@0.75 alongside Best RMSD, the distance from each experimental reference to its closest generated conformer, averaged across molecules with a measurable result. The hit rates describe how often an ensemble reaches a specified tolerance, whereas Best RMSD describes the distance of its closest approach without imposing a cutoff.

<!-- results-table-3:start -->
**Table 2. Recovery of the CASF bound conformation at the ChEMBL-count target.**

| Method | Hit@0.5 (%) ↑ | Hit@0.75 (%) ↑ | Best RMSD (Å) ↓ |
| --- | --- | --- | --- |
| Torsion raw | 35.1 | 66.0 | 0.679 |
| Torsion minimized | 46.8 | 68.1 | 0.644 |
| FlowR | 47.9 | 70.2 | 0.637 |
| RDKit minimized | 45.7 | 70.2 | 0.618 |
| RDKit raw | 43.6 | 66.0 | 0.601 |
| Torsional Diffusion | 56.4 | 73.4 | 0.582 |
| MCF drugs-L | 54.3 | 71.3 | 0.559 |
| NExT-Mol DMT-L | 53.2 | 76.6 | 0.543 |
| LoQI | **67.0** | 74.5 | 0.536 |
| Qwen 1.7B FSQ | 56.4 | **80.9** | **0.506** |
| ┄┄┄ | ┄┄┄ | ┄┄┄ | ┄┄┄ |
| ChEMBL3D-PB | <u>63.8</u> | <u>78.7</u> | <u>0.525</u> |

Hit@t is the percentage of all 94 entries with at least one retained conformer at RMSD ≤ t Å; missing or unmeasurable outputs count as failures. Best RMSD is the per-molecule minimum over retained conformers, averaged across entries with defined RMSD (92 for NExT-Mol, 93 for MCF, 94 for the others). It does not average distances over all generated conformers. Generators are sorted by decreasing unrounded Best RMSD, placing the lowest value immediately above the separate ChEMBL3D-PB reference. Bold and underlining identify the best and second-best distinct values, including ties and the stored reference. Candidate targets precede filtering. Source: [selected CASF records](results_tables/casf_selected.csv).
<!-- results-table-3:end -->

Qwen has the highest Hit@0.75 in the initial comparison, recovering 80.9% of core entries, followed by the stored ChEMBL3D ensemble at 78.7%. At the tighter 0.5 Å tolerance, LoQI has the highest hit rate at 67.0%, followed by ChEMBL3D at 63.8%; Qwen recovers 56.4%. Thus, reaching the experimental geometry within 0.75 Å does not necessarily translate into the most frequent matches within 0.5 Å. Qwen nevertheless has the lowest Best RMSD, although its difference from ChEMBL3D is modest, at 0.506 versus 0.525 Å.

Hit rate and Best RMSD need not rank methods in the same order. For example, minimized RDKit has a higher Hit@0.75 than raw RDKit, but a higher Best RMSD. Hit rate gives equal credit to all matches within the cutoff and does not measure how far the remaining entries miss it. Best RMSD retains the magnitude of those distances, including large misses, while averaging only over measurable entries. It still describes the closest member of each ensemble; neither metric establishes that all generated conformers lie near the experimental reference.

We then consider geometric diversity and proximity to the bound conformation together. Figure 1 places mean cluster count at 0.5 Å against Best RMSD, with more clusters and a lower RMSD defining the favorable region toward the lower right. Qwen combines the largest cluster count with the lowest Best RMSD in this comparison. The other generators with the highest cluster counts, random torsion sampling and FlowR, produce nearly as many clusters but struggle to approach the experimental references as closely: their Best RMSDs are higher. Broad geometric diversity therefore does not by itself ensure close recovery of the bound conformation. LoQI and ChEMBL3D approach the references with fewer clusters. Qwen's distinctive position therefore reflects its combination of broad geometric variation and close approaches to the observed structures, rather than a large RMSD advantage over every alternative. This combination makes it a candidate for dataset construction in the preliminary evaluation, subject to confirmation after benchmark molecules have been excluded from its training data.

<!-- results-figure-2:start -->
![Figure 1](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-2-diversity-recovery.png)

**Figure 1. Geometric diversity and proximity to the experimental bound conformation at the ChEMBL-count target.** Each point represents one evaluated pipeline or the stored ChEMBL3D-PB ensemble. The horizontal axis gives mean cluster count at 0.5 Å; the vertical axis gives Best RMSD, the per-molecule minimum over retained conformers averaged across entries with defined measurements. Dashed horizontal and vertical lines mark the ChEMBL3D-PB values. Points to the right and below these guides have more clusters and lower Best RMSD than the stored ensemble. Means use 92 entries for NExT-Mol, 93 for MCF, and 94 for the others. Axes show the observed region for readability. Retained ensemble sizes differ despite matched candidate targets. These are descriptive means without uncertainty intervals; Best RMSD does not describe every generated conformer.

[PDF](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-2-diversity-recovery.pdf) · [Editable SVG](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-2-diversity-recovery.svg)
<!-- results-figure-2:end -->

We next examine energy variation in the same retained ensembles. Figure 2 compares geometric cluster counts with the median of the per-molecule energy standard deviations, using a common MMFF94s calculation with unoptimized hydrogen coordinates and unchanged heavy-atom positions. The two panels describe the number of represented shapes and the spread of their force-field energies at the ChEMBL-count target.

<!-- results-figure-energy:start -->
![Figure 2](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-energy-dispersion.png)

**Figure 2. Geometric diversity and energy dispersion at the ChEMBL-count target.** (A) Mean geometric cluster count at 0.5 Å. (B) Median across molecules of the population standard deviation of conformer energies within each retained ensemble, in kcal/mol. Energies were recalculated with MMFF94s after rebuilding hydrogen coordinates without optimizing any atom. Each method uses the same measurable entries in both panels: 92 for NExT-Mol, 93 for MCF, and 94 for the others. The dashed line separates the stored ChEMBL3D-PB ensemble. Smaller energy SD indicates a narrower distribution; it does not establish lower absolute energy or the absence of high-energy conformers. Means, paired energy differences, and calculation completeness are reported in Supporting Table S9.

[PDF](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-energy-dispersion.pdf) · [Editable SVG](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-energy-dispersion.svg)
<!-- results-figure-energy:end -->

The methods with the largest cluster counts differ in energy spread. Qwen and FlowR have median energy standard deviations of 17.6 and 11.6 kcal/mol, respectively, compared with 40.5 kcal/mol for random torsion sampling and 26.5 kcal/mol for Torsional Diffusion. LoQI has a narrower energy distribution, with a median of 3.7 kcal/mol, while the minimized baselines have medians of 2.7–2.8 kcal/mol, alongside their lower cluster counts. The stored ChEMBL3D ensemble has a median energy standard deviation of 2.7 kcal/mol.

The mean across molecules is more affected by extreme energy spreads. For Torsional Diffusion, the mean per-molecule standard deviation is 1,628.4 kcal/mol, compared with its median of 26.5 kcal/mol. For Qwen, the corresponding summaries are 49.0 and 17.6 kcal/mol. The median therefore retains the differences in typical energy spread while reducing the influence of exceptional molecules. Supporting Table S9 reports both summaries and the corresponding energy levels. The median paired difference in mean energy relative to the same molecule's ChEMBL3D ensemble is +20.7 kcal/mol for Qwen and +0.8 kcal/mol for FlowR.

Matching the stored ChEMBL3D conformer counts establishes how the generators compare with an existing conformational resource. Constructing a larger dataset also requires understanding what additional sampling provides. We therefore extend the candidate target to 1,000 per molecule and examine whether the larger ensembles recover bound geometries missed at the smaller target, while increasing their geometric diversity.

All ten generation pipelines have higher Hit@0.75 and lower Best RMSD at the larger target. Qwen retains the highest hit rate and lowest Best RMSD, increasing recovery from 80.9% to 91.5% and reducing Best RMSD from 0.506 to 0.307 Å. Its mean number of clusters at 0.5 Å increases from 55.3 to 284.5. FlowR and random torsion sampling produce the largest cluster counts at 1,000 candidates, but their recovery remains below that of Qwen. The increase in geometric diversity therefore does not translate into the same recovery across generators.

The stored ChEMBL3D ensemble remains unchanged, with a hit rate of 78.7%. Additional sampling extends recovery beyond the available stored conformers, but the comparison does not describe what ChEMBL3D would achieve with an equally enlarged ensemble. The two candidate targets also yield different retained counts across methods and do not represent equal computational costs.

<!-- results-table-4:start -->
**Table 3. Recovery and geometric diversity at the two candidate targets.**

| Method | 1,000: Best RMSD (Å) ↓ | 1,000: Hit@0.75 (%) ↑ | 1,000: clusters at 0.5 Å ↑ | ChEMBL-count: Best RMSD (Å) ↓ | ChEMBL-count: Hit@0.75 (%) ↑ | ChEMBL-count: clusters at 0.5 Å ↑ |
| --- | --- | --- | --- | --- | --- | --- |
| Torsion minimized | 0.535 | 79.8 | 60.0 | 0.644 | 68.1 | 21.9 |
| RDKit minimized | 0.512 | 78.7 | 52.1 | 0.618 | 70.2 | 21.7 |
| Torsion raw | 0.477 | 78.7 | <u>296.8</u> | 0.679 | 66.0 | <u>54.2</u> |
| RDKit raw | 0.437 | 80.9 | 114.6 | 0.601 | 66.0 | 34.3 |
| FlowR | 0.414 | 86.2 | **298.1** | 0.637 | 70.2 | 54.1 |
| NExT-Mol DMT-L | 0.386 | 86.2 | 106.5 | 0.543 | 76.6 | 31.6 |
| Torsional Diffusion | 0.371 | <u>89.4</u> | 187.7 | 0.582 | 73.4 | 40.9 |
| MCF drugs-L | 0.366 | 86.2 | 110.3 | 0.559 | 71.3 | 31.2 |
| LoQI | <u>0.339</u> | 86.2 | 84.5 | 0.536 | 74.5 | 25.0 |
| Qwen 1.7B FSQ | **0.307** | **91.5** | 284.5 | **0.506** | **80.9** | **55.3** |
| --- | --- | --- | --- | --- | --- | --- |
| ChEMBL3D-PB | — | — | — | <u>0.525</u> | <u>78.7</u> | 21.7 |

Hit@0.75 uses all 94 entries, with missing outputs counted as failures. Best RMSD is the mean of the per-entry minimum RMSD; cluster counts are also averaged over entries with defined measurements. At the 1,000-candidate target, these means use 93 entries for MCF and NExT-Mol and 94 for the others; at the ChEMBL-count target, they use 93 for MCF, 92 for NExT-Mol, and 94 for the others. Generators are sorted by decreasing unrounded Best RMSD at 1,000 candidates. Bold and underlining identify the best and second-best distinct values in each column, including ties and the stored reference where reported. ChEMBL3D-PB is the same stored ensemble and has no 1,000-candidate result. Candidate targets precede filtering and do not match retained counts or computational cost.

Source: [budget comparison](results_tables/sampling_budget_comparison.csv); [per-entry records](results_tables/sampling_budget_per_entry.csv).
<!-- results-table-4:end -->

<!-- results-figure-3:start -->
![Figure 3](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-3-sampling-budget.png)

**Figure 3. Recovery and diversity at the two candidate targets.** Open circles denote the ChEMBL-count target and filled circles the 1,000-candidate target. Each row follows the same method between the two targets. The stored ChEMBL3D-PB ensemble is shown once as a hexagon. Recovery at 0.75 Å uses all 94 CASF entries; mean cluster counts at 0.5 Å use defined measurements. Rows follow Table 3, ordered by decreasing Best RMSD at 1,000 candidates; the dashed line separates the stored reference. Lines connect observed endpoints and do not represent random-subsampling curves or intermediate measurements. Targets precede PoseBusters filtering and do not imply equal retained ensemble sizes or computational costs.

[PDF](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-3-sampling-budget.pdf) · [Editable SVG](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-3-sampling-budget.svg)
<!-- results-figure-3:end -->

Recovery establishes whether an ensemble contains at least one close match. The fraction of retained conformers within 0.75 Å of the reference addresses a complementary question: how frequently the method samples near that observed geometry. This distinction matters for a training resource, where the distribution of conformers contributes to the data presented to a model.

[If included, report this conformer-level fraction separately from ligand-level Hit@0.75, with its denominator and aggregation rule. Insert its interpretation only after calculation on the final valid pools.]

Overall recovery can conceal substantial differences between molecules. We therefore compare recovery across groups defined by rotatable-bond and heavy-atom counts to assess how aggregate performance changes with molecular flexibility and size.

<!-- results-figure-4:start -->
![Figure 4](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-4-size-flexibility.png)

**Figure 4. Recovery by molecular size and flexibility.** Generated methods use the 1,000-candidate target; ChEMBL3D-PB uses its stored ensemble. Every cell gives recovered/total CASF entries, with the color indicating the corresponding percentage. Dashed outlines mark groups with fewer than five entries. Heavy-atom and rotatable-bond groups are analyzed separately and do not isolate independent effects of size and flexibility. Sparse groups cannot establish stable method rankings.

[PDF](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-4-size-flexibility.pdf) · [Editable SVG](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-4-size-flexibility.svg)
<!-- results-figure-4:end -->

The core panel contains relatively few molecules in the largest and most flexible groups. Its stratified results therefore identify observed successes and failures but cannot establish a stable ranking in every regime. These differences are relevant to dataset construction because a single average recovery rate can hide poor coverage of precisely the molecules that require the most extensive conformational sampling.

The CASF analysis asks whether an ensemble contains a close match to one observed conformation for each entry. A broader training resource should also represent variation among the conformations that the same molecule adopts in different experimental structures. We therefore extend the evaluation to a panel of 23 drug-like molecules with multiple observed bound conformations selected from PLINDER release 2024-06/v2 by physicochemical screening, PDB-frequency ranking, and manual review (Methods; Appendix A). For each molecule, the same protein-independent ensemble is compared with all of its experimental references. This tests coverage of several observed geometries without generating a separate ensemble for each complex. [PLINDER](https://www.biorxiv.org/content/10.1101/2024.07.17.603955v3)

Reference coverage and sample precision distinguish how many experimental observations an ensemble covers from how frequently it samples near them. Table 4 reports these two measurements together with the matching distances defined in Methods.

<!-- results-table-5:start -->
**Table 4. Exploratory coverage and precision for the 23-molecule PLINDER panel.**

| Method | Mean supplied N | COV-R (%) | COV-P (%) | MAT-R (Å) | MAT-P (Å) | PB pass (%) |
| --- | --- | --- | --- | --- | --- | --- |
| LoQI | 1000.0 | 91.2 | 59.0 | 0.421 | 0.816 | 90.9 |
| Torsional Diffusion | 1000.0 | 78.6 | 42.6 | 0.455 | 1.189 | 83.4 |
| MCF drugs-L | 1000.0 | 82.7 | 39.1 | 0.481 | 1.215 | 63.7 |
| NExT-Mol DMT-L | 1000.0 | 82.7 | 36.9 | 0.464 | 1.172 | 62.2 |
| FlowR | 961.7 | 81.3 | 36.6 | 0.522 | 1.306 | 78.5 |
| Qwen 1.7B FSQ | 999.2 | 91.0 | 46.8 | 0.375 | 1.000 | 83.7 |

COV-R and COV-P use a strict RMSD < 0.75 Å cutoff in the existing evaluation; MAT-R and MAT-P are continuous distances. Coverage, precision, matching distances, and supplied counts are means over 23 molecules, with 2,450 reference records in total. PB pass is the separately measured pooled pass fraction. The coverage and distance columns use supplied pools without a common PB filter, and the evaluators differ in RMSD-failure handling. They therefore remain descriptive, pending harmonized evaluation; unmatched conformers are not necessarily physically inaccessible.

Sources: [molecule-level records](results_tables/druglike_selected.csv) and [unrounded table values](results_tables/druglike_means.csv).


<!-- results-table-5:end -->

<!-- results-figure-5:start -->
![Figure 5](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-5-multiple-references.png)

**Figure 5. Coverage of multiple bound references and precision relative to those observations.** (A) Molecule-averaged COV-R and COV-P for six learned generators, using the existing strict RMSD < 0.75 Å criterion. Axes are restricted to the observed range for readability. (B) Paired reference coverage for Qwen and LoQI across all 23 molecules; the diagonal denotes equal coverage. Coincident points are grouped and labeled, including 11 molecules with complete coverage by both methods. Imatinib and actinonin are the examples retained from the earlier draft. These are descriptive supplied-pool results with differing RMSD-failure conventions, not a comparison after common validity filtering. No uncertainty intervals are shown; paired uncertainty is documented in the archived findings.

[PDF](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-5-multiple-references.pdf) · [Editable SVG](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-5-multiple-references.svg)
<!-- results-figure-5:end -->

In the existing evaluation, LoQI and Qwen achieve similar mean reference coverage, while LoQI has higher observed-reference precision. These preliminary comparisons describe the supplied pools; common validity filtering and consistent treatment of failed RMSD calculations remain necessary for the final evaluation. Thus, similar coverage can be obtained from ensembles that differ in the proportion of samples lying near the available experimental structures. This distinction qualifies the single-reference recovery result and is relevant when considering how generated conformers will populate a training dataset. Experimental references remain incomplete, however, so an unmatched conformer cannot automatically be classified as physically inaccessible or irrelevant to binding.

[Recompute the drug-panel comparison using a common validity filter and RMSD-failure policy. The PLINDER release and molecule-selection criteria are documented in Methods and Appendix A; complete the description of reference-coordinate extraction and handling of repeated or nearly identical structures.]

Together, experimental coverage, geometric variation, and structural and energetic plausibility inform the selection of a generation procedure for the larger dataset. Our current plan is to apply Qwen to molecules drawn from ChEMBL3D, subject to confirmation by the cleaned evaluation. The dataset will be characterized in terms of molecular coverage, conformers per molecule, validity yield, geometric diversity, and energy distributions. These measurements will establish what the selected procedure produces at scale and whether the resulting resource retains the properties that motivated its selection. The intended uses span structure-based drug design and learning molecular geometry, including pharmacophore-conditioned generation; performance in downstream applications remains to be evaluated separately.

<!-- results-table-6:start -->
**Table 5. Dataset-release reporting status (author working table).**

| Release item | Current status |
| --- | --- |
| Molecular source | ChEMBL3D molecules are planned; final selection pending |
| Generator | Qwen is the planned candidate, subject to the evaluation after benchmark exclusion |
| Molecule and conformer counts | Pending generation and characterization |
| Conformers per molecule and validity yield | Pending; benchmark counts above are not release counts |
| Molecular size and flexibility distributions | Pending final molecular selection |
| Geometric diversity and energy distributions | Pending analysis of the released pools |
| Generation settings and checkpoint | Final configuration pending |
| Public access and version identifier | Pending release |

This table records unfinished release items; it is not a measured dataset result and should be replaced after generation and characterization.
<!-- results-table-6:end -->

# Discussion

The central purpose of this study is to identify a practical procedure for constructing a large conformational dataset with coverage of experimentally observed ligand geometries. The benchmark separates several properties that are often combined under the general description of conformer quality. Structural validity determines whether a sample satisfies the specified identity and geometry requirements. Geometric diversity describes the variation represented in an ensemble. Energy analysis characterizes that variation under a common molecular model, while experimental-conformation recovery assesses whether the ensemble includes shapes observed in protein complexes.

The preliminary results support the value of considering these properties together. Broad geometric exploration does not guarantee recovery, and strong recovery by the closest member of an ensemble does not establish that most samples lie near observed conformations. The comparison between single-reference recovery and multi-reference coverage makes this distinction explicit. It also explains why the preferred generator may depend on whether the immediate objective is to recover at least one useful conformation, cover several experimental conformations, or concentrate sampling near those observations.

These evaluations inform dataset construction without establishing that every generated conformer is biologically relevant. The experimental collections sample only part of the accessible conformational space, and the force-field calculations do not represent binding free energies. The results are also conditional on the selected molecules, generation settings, and sampling budgets. Conclusions about Qwen’s generalization will depend on the evaluation after benchmark molecules have been excluded from its training data.

The proposed public release is intended to provide a reproducible resource for structure-based drug design and learning molecular geometry, with pharmacophore-conditioned generation among its potential applications. Its demonstrated properties will be the structural coverage, diversity, and validity measured in this study. Whether training on this resource improves a particular downstream model remains a separate question for future work.
