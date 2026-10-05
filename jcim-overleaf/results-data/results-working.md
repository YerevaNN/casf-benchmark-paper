# Results

> **Data status:** Tables use the corrected findings archived on 24 September 2026. Qwen values belong to the initial 1.7B FSQ checkpoint at step 47,023; they do not describe the pending evaluation after benchmark exclusions. Energy results and release counts remain pending.

A useful conformational resource should represent a range of molecular shapes while retaining coverage of geometries observed in protein complexes. We begin by examining the variation present in the generated ensembles, before asking whether that variation includes the experimental bound conformations. The first comparison uses the 94-entry CASF core panel at the ChEMBL-count target, placing each generator alongside the stored ChEMBL3D ensemble for the same molecules. Table 2 summarizes geometric diversity through cluster counts at 1.0 Å, together with the retained ensemble sizes.

<!-- results-table-2:start -->
**Table 2. Geometric diversity at the ChEMBL-count target on the 94-entry core panel.**

| Method | Mean retained conformers | Mean clusters at 1.0 Å |
| --- | --- | --- |
| ChEMBL3D-PB | 81.2 | 7.7 |
| RDKit raw | 81.2 | 14.7 |
| RDKit minimized | 81.2 | 11.0 |
| Torsion raw | 71.1 | 21.8 |
| Torsion minimized | 81.2 | 11.7 |
| LoQI | 81.2 | 9.2 |
| Torsional Diffusion | 69.3 | 17.5 |
| MCF drugs-L | 66.1 | 14.2 |
| NExT-Mol DMT-L | 66.4 | 14.2 |
| FlowR | 72.7 | 29.1 |
| Qwen 1.7B FSQ | 78.3 | 22.7 |

ChEMBL3D-PB is the available stored ensemble. Retained counts are averaged over all 94 entries, including empty outputs. Cluster counts are means over entries with defined clustering measurements; these are geometric clusters, not energy basins. Source: [selected CASF records](results_tables/casf_selected.csv).
<!-- results-table-2:end -->

The initial comparison shows that methods differ substantially in the geometric variation they produce. Broad torsional perturbation generates many distinct shapes, while learned generators also vary considerably in how widely they explore conformational space. FlowR and Qwen illustrate that substantial diversity can arise from different generation procedures, whereas LoQI provides a contrasting ensemble shaped by its low-energy objective. These observations describe the geometry of the sampled pools; they do not yet establish whether the additional shapes are energetically plausible or close to experimentally observed conformations.

[Update the method-specific observations and their magnitudes from the final ChEMBL-count evaluation; do not transfer rankings from the 1,000-candidate comparison to this target.]

We next examine the energies associated with this geometric variation. The planned MMFF94s analysis will assess whether broad geometric sampling is accompanied by a substantial high-energy tail, using within-molecule comparisons and the robust summaries described in Methods.

Energy dispersion provides information that cluster counts alone cannot supply. An ensemble may contain many geometrically different conformers while assigning substantial sampling effort to highly strained structures. Conversely, a narrow energy distribution does not necessarily imply limited geometric variation, since distinct conformations can have similar energies. We therefore interpret this analysis as a measure of energetic plausibility and dispersion, rather than as a count of occupied energy basins. The intended dataset should preserve meaningful structural variation without being dominated by an extreme high-energy tail.

[Here goes the updated MMFF94s energy figure, showing within-molecule energy comparisons, robust summaries of energy spread, and the high-energy tail. Recompute this analysis using the corrected molecular identities and final generated pools before inserting a result-specific interpretation.]

Having characterized the ensembles themselves, we assess recovery of the experimental bound geometries at the ChEMBL-count target. Table 3 reports ligand-level Hit@0.75 and the mean minimum RMSD, together with retained ensemble sizes and the number of measurable entries.

<!-- results-table-3:start -->
**Table 3. Recovery of the CASF bound conformation at the ChEMBL-count target.**

| Method | Recovered entries | Hit@0.75 (%) | Mean minimum RMSD (Å) | Entries with measurable RMSD | Mean retained conformers |
| --- | --- | --- | --- | --- | --- |
| ChEMBL3D-PB | 74/94 | 78.7 | 0.525 | 94 | 81.2 |
| RDKit raw | 62/94 | 66.0 | 0.601 | 94 | 81.2 |
| RDKit minimized | 66/94 | 70.2 | 0.618 | 94 | 81.2 |
| Torsion raw | 62/94 | 66.0 | 0.679 | 94 | 71.1 |
| Torsion minimized | 64/94 | 68.1 | 0.644 | 94 | 81.2 |
| LoQI | 70/94 | 74.5 | 0.536 | 94 | 81.2 |
| Torsional Diffusion | 69/94 | 73.4 | 0.582 | 94 | 69.3 |
| MCF drugs-L | 67/94 | 71.3 | 0.559 | 93 | 66.1 |
| NExT-Mol DMT-L | 72/94 | 76.6 | 0.543 | 92 | 66.4 |
| FlowR | 66/94 | 70.2 | 0.637 | 94 | 72.7 |
| Qwen 1.7B FSQ | 76/94 | 80.9 | 0.506 | 94 | 78.3 |

Recovery requires at least one retained conformer at RMSD ≤ 0.75 Å. Missing or unmeasurable outputs count as failures among all 94 entries. Mean minimum RMSD uses only entries with defined RMSD; its denominator is shown explicitly. Candidate targets precede filtering. Source: [selected CASF records](results_tables/casf_selected.csv).
<!-- results-table-3:end -->

The preliminary results show that the stored ChEMBL3D ensembles already recover a substantial fraction of the experimental conformations. At the ChEMBL-count target, the selected Qwen model performs similarly to this baseline on the core panel. The existing dataset therefore provides a meaningful reference for assessing whether generated ensembles offer useful experimental coverage at a comparable candidate count. This comparison does not by itself establish an advantage for Qwen over the stored ensemble.

We then consider diversity and recovery together. The initial results show that a larger number of geometric clusters does not consistently correspond to better recovery of bound conformations. Some methods explore a broad range of shapes but still miss experimental references that are recovered by less diverse ensembles. The position of each method in the joint comparison indicates whether its geometric variation is accompanied by recovery of the experimental references. Qwen is a candidate for dataset construction because it combines geometric variation with recovery in the preliminary evaluation, but this position must be confirmed after exclusion of benchmark molecules from its training data. The relevant distinction is therefore not only how much conformational variation a method produces, but whether that variation includes the regions represented by experimental observations.

<!-- results-figure-2:start -->
![Figure 2](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-2-diversity-recovery.png)

**Figure 2. Geometric diversity and experimental recovery at the ChEMBL-count target.** Each point represents one evaluated pipeline or the stored ChEMBL3D-PB ensemble. Clustering uses a 1.0 Å radius, whereas recovery requires at least one retained conformer at RMSD ≤ 0.75 Å. Cluster counts are averaged over defined measurements; recovery uses all 94 entries, including failures. Retained ensemble sizes differ despite matched candidate targets. These are descriptive method summaries, not a causal analysis of diversity.

[PDF](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-2-diversity-recovery.pdf) · [Editable SVG](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-2-diversity-recovery.svg)
<!-- results-figure-2:end -->

Constructing a large conformational resource also requires understanding what additional sampling provides. We therefore extend the comparison to the target of 1,000 candidates per molecule. This tests whether the relationship between diversity and recovery persists as the ensembles grow and whether additional samples reach experimental geometries missed at the smaller target.

Increasing the candidate allowance improves recovery in the initial evaluation. Qwen recovers additional bound conformations at the 1,000-candidate target, including conformations absent from the stored ChEMBL3D pools. This suggests that additional generated sampling can extend the structural coverage available for the same molecules. The comparison with ChEMBL3D concerns its available stored conformers: its ensemble does not grow when the generated methods receive a larger sampling allowance. The high-budget comparison therefore measures the additional coverage supplied by generation rather than superiority over an equally extended ChEMBL3D sampling procedure. Differences in retained counts also mean that the candidate targets should not be interpreted as identical valid ensemble sizes or computational costs.

<!-- results-table-4:start -->
**Table 4A. Recovery at the two candidate targets.**

| Method | ChEMBL-count: Hit (%) | ChEMBL-count: retained N | ChEMBL-count: min. RMSD (Å) | 1,000: Hit (%) | 1,000: retained N | 1,000: min. RMSD (Å) |
| --- | --- | --- | --- | --- | --- | --- |
| ChEMBL3D-PB | 78.7 | 81.2 | 0.525 | — | — | — |
| RDKit raw | 66.0 | 81.2 | 0.601 | 80.9 | 999.9 | 0.437 |
| RDKit minimized | 70.2 | 81.2 | 0.618 | 78.7 | 1000.0 | 0.512 |
| Torsion raw | 66.0 | 71.1 | 0.679 | 78.7 | 896.5 | 0.477 |
| Torsion minimized | 68.1 | 81.2 | 0.644 | 79.8 | 1000.0 | 0.535 |
| LoQI | 74.5 | 81.2 | 0.536 | 86.2 | 999.3 | 0.339 |
| Torsional Diffusion | 73.4 | 69.3 | 0.582 | 89.4 | 881.4 | 0.371 |
| MCF drugs-L | 71.3 | 66.1 | 0.559 | 86.2 | 866.8 | 0.366 |
| NExT-Mol DMT-L | 76.6 | 66.4 | 0.543 | 86.2 | 871.9 | 0.386 |
| FlowR | 70.2 | 72.7 | 0.637 | 86.2 | 920.0 | 0.414 |
| Qwen 1.7B FSQ | 80.9 | 78.3 | 0.506 | 91.5 | 972.6 | 0.307 |

Hit uses the 0.75 Å cutoff and all 94 entries. Retained N and minimum RMSD are means, with the same denominator conventions as Table 3. At the 1,000-candidate target, minimum RMSD is defined for 93 entries for MCF and NExT-Mol and all 94 for the other pipelines. ChEMBL3D-PB is shown once because its stored ensemble does not grow. These comparisons do not match retained counts or computational cost.

**Table 4B. Geometric diversity in the larger ensembles.**

| Method | Mean retained conformers | Mean clusters at 1.0 Å |
| --- | --- | --- |
| ChEMBL3D-PB | 81.2 | 7.7 |
| RDKit raw | 999.9 | 35.9 |
| RDKit minimized | 1000.0 | 20.2 |
| Torsion raw | 896.5 | 74.2 |
| Torsion minimized | 1000.0 | 25.1 |
| LoQI | 999.3 | 20.2 |
| Torsional Diffusion | 881.4 | 54.9 |
| MCF drugs-L | 866.8 | 34.0 |
| NExT-Mol DMT-L | 871.9 | 31.5 |
| FlowR | 920.0 | 114.4 |
| Qwen 1.7B FSQ | 972.6 | 72.8 |

Generated methods use the 1,000-candidate target; ChEMBL3D-PB remains the same stored ensemble. Cluster means use defined measurements, as in Table 2. Source for both panels: [selected CASF records](results_tables/casf_selected.csv).


<!-- results-table-4:end -->

<!-- results-figure-3:start -->
![Figure 3](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-3-sampling-budget.png)

**Figure 3. Recovery and diversity at the two candidate targets.** Open circles denote the ChEMBL-count target and filled circles the 1,000-candidate target. Each row follows the same method between the two targets. The stored ChEMBL3D-PB ensemble is shown once as a hexagon. Recovery uses all 94 CASF entries; cluster means use defined measurements. Lines connect observed endpoints and do not represent random-subsampling curves or intermediate measurements. Targets precede PoseBusters filtering; retained counts are given in Table 4.

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

Reference coverage and sample precision distinguish how many experimental observations an ensemble covers from how frequently it samples near them. Table 5 reports these two measurements together with the matching distances defined in Methods.

<!-- results-table-5:start -->
**Table 5. Exploratory coverage and precision for the 23-molecule PLINDER panel.**

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

Together, experimental coverage, geometric variation, and structural and energetic plausibility inform the selection of a generation procedure for the larger dataset. Our current plan is to apply Qwen to molecules drawn from ChEMBL3D, subject to confirmation by the cleaned evaluation. The dataset will be characterized in terms of molecular coverage, conformers per molecule, validity yield, geometric diversity, and energy distributions. These measurements will establish what the selected procedure produces at scale and whether the resulting resource retains the properties that motivated its selection. The principal intended use is pharmacophore-conditioned learning, with other applications in molecular geometry and representation learning; downstream training performance remains to be evaluated separately.

<!-- results-table-6:start -->
**Table 6. Dataset-release reporting status (author working table).**

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

The proposed public release is intended to provide a reproducible resource for learning molecular geometry, particularly for pharmacophore-conditioned generation. Its demonstrated properties will be the structural coverage, diversity, and validity measured in this study. Whether training on this resource improves a particular downstream model remains a separate question for future work.
