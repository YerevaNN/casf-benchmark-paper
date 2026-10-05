# Methods

Working draft, 5 October 2026. Table and figure numbers are shared with the Results draft. Bracketed author notes identify unfinished evidence and are not manuscript prose. This draft describes the archived evaluation and distinguishes it from the pending Qwen, energy, and dataset-release work.

## Evaluation design

We evaluate how well finite, protein-independent conformer ensembles cover experimentally observed bound geometries, with the practical aim of selecting a generation procedure for a large conformational dataset. Each method generates an ensemble for a molecule without receiving a protein or binding-site structure. We examine geometric variation and energy alongside experimental recovery, before extending the comparison from one reference per CASF entry to multiple bound structures of the same molecule. Recovery measures similarity of the internal ligand geometry after alignment; it does not establish placement within a binding pocket or recovery of protein–ligand interactions.

## Experimental references and the ChEMBL3D comparison

We compare generated and existing conformer ensembles using CASF-2016, a curated dataset of 285 protein–ligand crystal complexes. The bound ligand geometries provide experimental references for assessing whether a protein-independent ensemble includes an observed conformation. [CASF-2016](https://doi.org/10.1021/acs.jcim.8b00545)

To compare generation with an existing dataset, we identify the CASF ligands that are also represented in ChEMBL3D. CASF ligand MOL2 files are parsed with RDKit and matched to the ChEMBL3D index by canonical heavy-atom isomeric SMILES. Matching therefore preserves molecular connectivity and stereochemistry. Because several stereoisomers can share a ChEMBL3D parent identifier, we resolve the matching topology SDF record using stereochemistry reconstructed from its coordinates, comparing both isomeric SMILES and the tetrahedral and double-bond stereochemistry layers of InChI. Stored conformers and their counts are restricted to the matching stereoisomer. We do not choose an isomer by its RMSD to the experimental reference.

Entries are retained when the ligand can be parsed, a matching topology is available, and the topology contains at least one eligible torsion. Torsion eligibility uses the SMARTS pattern `[!$(*#*)&!D1]-!@[!$(*#*)&!D1]`, which differs from the rotatable-bond descriptor used for the flexibility analysis. These criteria produce a core panel of 94 complex entries representing 94 ChEMBL parent compounds. Applying the corresponding procedure to the larger reference collection produces 1,236 entries representing 1,042 parent compounds. The collections share 17 parent compounds but no complex identifiers. Core provides the main comparison; the larger collection is reported in the Supporting Information.

For each molecule, the stored ChEMBL3D ensemble establishes how much of the experimental geometry is already represented in an existing computed dataset. This comparison distinguishes coverage already available in the reference resource from coverage that can be added by a generation procedure. Because the comparison is restricted to molecules represented in ChEMBL3D, its conclusions concern conformational coverage within this intersection rather than coverage of the entire CASF chemical space. The analysis permits a deterministic cap of 2,000 stored conformers for the designated large reference entry, 1tlp, when its available ensemble exceeds that size.

The complementary experiment uses an exploratory panel of 23 ligand identities selected from PLINDER release 2024-06/v2, an annotated resource of protein–ligand complexes from the Protein Data Bank. Records in the ligand-per-system table were grouped by full InChIKey so that different stereoisomers remained separate. Of 51,280 starting identities, 25,392 passed an RDKit-based physicochemical screen requiring molecular weight 150–650 Da, cLogP −1 to 6, and limits on polarity, hydrogen bonding, flexibility, charge, and elemental composition. The full criteria are given in Appendix Table A1. [PLINDER](https://doi.org/10.1101/2024.07.17.603955)

Survivors were ranked by the number of distinct PDB entries, with PLINDER-system count as the secondary ordering. Distinct PDB entries were preferred because one multimeric structure can contribute many systems. Manual review then selected recognizable drugs and medicinal-chemistry-like protein ligands, excluding common crystallization additives, metabolites, cofactors, and other compounds outside the intended scope. The resulting Plinder-23 panel is a manually curated exploratory set, not an official PLINDER split or a purely algorithmic selection. Exact reproduction uses the saved identities. The molecule list and panel statistics are provided in [Appendix A](appendix_working.md).

The evaluated panel contains 2,450 reference conformers, with 23–716 per molecule. Each generated ensemble is compared with all supplied references for that molecule. The saved PDB and system frequencies describe selection, whereas these conformer counts describe the evaluated coordinate collection. Coverage concerns reference records rather than a separately defined set of distinct conformational states.

[Author note: confirm the source release and construction of the larger CASF reference collection. For Plinder-23, the release and identity-selection procedure are now documented; reference-coordinate extraction and preprocessing still need to be specified, including protonation and tautomer handling and treatment of repeated or nearly identical experimental structures.]

<!-- results-figure-1:start -->
![Figure 1](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-1-datasets.png)

**Figure 1. Experimental evaluation panels.** CASF-2016 is matched to ChEMBL3D by molecular identity and stereochemistry, followed by eligibility filtering, to obtain 94 core entries. The stored ChEMBL3D ensemble is the computed comparison resource. The separate 23-molecule PLINDER panel assesses several references per molecule. The larger 1,236-entry collection remains in Supporting Information. Arrows describe evaluation design, not a quantified exclusion funnel. For PLINDER 2024-06/v2, physicochemical screening reduced 51,280 identities to 25,392; ranking by distinct PDB entries and manual review yielded the 23-molecule panel (Methods; Appendix A). The supplied-pool comparison does not yet use the common CASF validity filter; reference-coordinate preprocessing remains to be documented.

[PDF](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-1-datasets.pdf) · [Editable SVG](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-1-datasets.svg)
<!-- results-figure-1:end -->

## Conformer generators

The generators represent different approaches to conformational sampling. The main comparison contains four classical pipelines, five external learned generators, and one Qwen model. Table 1 summarizes their representations, training provenance, initialization, and processing. These methods differ in training data and inference procedures as well as architecture, so the comparison evaluates the resulting generation pipelines rather than isolating any single modeling choice.

RDKit ETKDGv3 provides an established distance-geometry baseline incorporating experimental torsional preferences. Conformers are embedded from a coordinate-free copy of the matched molecular topology, using random starting coordinates, enforced chirality, and no RMSD pruning. The raw pipeline rejects nonfinite coordinates without a separate pre-generation clash filter. The minimized pipeline applies MMFF94s for up to 500 iterations. Optimization status 0 or 1 is accepted, so reaching the iteration limit does not require rejection. Raw and minimized variants use independently sampled pools; they are not a paired before-and-after minimization experiment. [ETKDG](https://doi.org/10.1021/acs.jcim.0c00025)

The random torsion-perturbation baseline explores broad torsional variation with limited distributional guidance. It starts from the matched ChEMBL3D topology's stored geometry and independently perturbs every eligible torsion by a uniform displacement between −120° and +120°. Other aspects of the initial geometry are inherited. Candidates are rejected when a nonbonded heavy-atom distance is below 0.7 times the distance-geometry lower bound; bonded and 1,3 pairs are excluded from this check. The minimized variant adds MMFF94s minimization and repeats the clash check. These pipelines accumulate candidates toward the requested pool size before the common validity filter is applied.

[Author note: the shared reference loader can fall back to a stereochemically matching CASF MOL2 when the ChEMBL3D topology cannot be loaded. This fallback can also supply the torsion seed, not only the validity reference. Its use in the archived runs has not been counted. Audit it before stating that all torsion seeds were independent of the experimental geometry; rerun affected entries from an independent seed if necessary.]

The learned methods provide complementary sampling strategies. Torsional Diffusion learns a distribution over torsion angles, while LoQI targets low-energy molecular conformations using a model trained on ChEMBL3D. MCF drugs-L and NExT-Mol DMT-L provide additional coordinate-generation approaches. We evaluate their saved outputs through the benchmark adapters. FlowR uses the ligand-only FLOWR.root v2.2 checkpoint with the input atom and bond graph fixed during coordinate generation. Its inspected adapter uses harmonic graph inpainting, stereochemistry checks, and up to 200 force-field cleanup iterations with MMFF and a UFF fallback. A reflected structure may be retained if reflection restores the requested stereoisomer; remaining mismatches are rejected. Thus, an output labeled raw can already include generator-specific processing. [Torsional Diffusion](https://arxiv.org/abs/2206.01729), [LoQI](https://doi.org/10.26434/chemrxiv-2025-k4h7v)

Qwen represents autoregressive generation of a molecular coordinate representation conditioned on molecular SMILES. The archived comparison uses the 1.7-billion-parameter FSQ model at pretraining step 47,023, identified as `qwen_1p7b_fsq_bigdata_step47023`, for both CASF and the 23-molecule experiment. Its initial training corpus was broadly assembled, and these results do not establish generalization after exclusion of benchmark molecules. The planned evaluation will use a model trained after excluding CASF and PLINDER benchmark molecules. External pretrained checkpoints are evaluated as supplied, with their documented training provenance and any overlap that can be checked.

<!-- TODO: insert the verified Qwen coordinate-tokenizer citation. -->
[TODO: add the Qwen coordinate-tokenizer citation, final training description and exclusion rules. Record checkpoint selection, hashes, dependency versions, decoding settings, seeds, retries, and executed external inference configurations. Current launch defaults do not by themselves establish the settings used for historical outputs.]

<!-- results-table-1:start -->
**Table 1. Evaluated generation pipelines and stored comparison ensemble.**

| Method | Representation / sampling | Training data and checkpoint | Input / initialization | Processing before common CASF filtering |
| --- | --- | --- | --- | --- |
| ChEMBL3D-PB | Stored optimized coordinates | Reference dataset | Stored ensemble for the matched stereoisomer | PoseBusters filter; no new generation |
| RDKit raw | ETKDGv3 distance geometry | No learned checkpoint | Molecular graph; random starting coordinates | Finite-coordinate check; no added MMFF94s minimization |
| RDKit minimized | ETKDGv3 distance geometry | No learned checkpoint | Same embedding procedure; independently sampled pool | MMFF94s, up to 500 iterations |
| Torsion raw | Uniform torsion perturbations | No learned checkpoint | Matched ChEMBL3D seed geometry; ±120° per eligible torsion | Steric clash filtering |
| Torsion minimized | Uniform torsion perturbations | No learned checkpoint | Same seed and perturbation procedure; independently sampled pool | Clash filtering, MMFF94s, repeat clash check |
| LoQI | Equivariant coordinate diffusion | ChEMBL3D; upstream checkpoint, exact hash pending | Molecular graph; exact initialization to document | Imported pipeline; postprocessing audit pending |
| Torsional Diffusion | Diffusion over rotatable torsions | GEOM-DRUGS checkpoint; exact hash pending | RDKit initial geometry | Imported pipeline; postprocessing audit pending |
| MCF drugs-L | Coordinate diffusion with molecular conformer fields | GEOM-DRUGS; drugs-L checkpoint, exact hash pending | Molecular graph; exact initialization to document | Imported pipeline; postprocessing audit pending |
| NExT-Mol DMT-L | Molecular language model with a coordinate diffusion head | GEOM-DRUGS conformer training; ZINC-15 language-model pretraining; DMT-L, exact hash pending | Molecular representation; exact initialization to document | Imported pipeline; postprocessing audit pending |
| FlowR | Coordinate flow matching with a fixed input graph | FLOWR.root v2.2; exact training composition and hash to confirm | Ligand-only harmonic graph inpainting | Up to 200 cleanup steps (MMFF, UFF fallback); stereochemistry checks |
| Qwen 1.7B FSQ | Autoregressive coordinate-token generation | Broadly assembled corpus; exact composition pending; pretraining step 47,023 | Molecular SMILES; FSQ representation | Coordinate decoding; exact inference settings to document |

All CASF ensembles subsequently undergo the specified PoseBusters checks. “Raw” does not imply an absence of generator-specific processing. Metadata describe the archived evaluation; settings and hashes marked pending need confirmation from the executed runs. The Qwen results below use the initial step-47,023 checkpoint, not a model retrained after benchmark exclusions. Its coordinate-tokenizer citation remains TODO. Training-set overlap has not been ruled out for the external checkpoints.

Sources: [generation procedures](generation_methods.md), [model catalogue](generator_models_catalog.md), [run registry](../src/casf_benchmark/config/generation_runs.yaml), and [manuscript inference audit](/mnt/weka/mbedrosian/code/casf-benchmark-paper/jcim-overleaf/METHODS_NOTES.md). Older catalogue descriptions of Qwen training are not used to infer its final training recipe.
<!-- results-table-1:end -->

## Sampling targets and retained ensembles

We use two sampling targets to reflect different aspects of dataset construction. At the ChEMBL-count target, each method is assigned the stored conformer count for the corresponding ChEMBL3D stereoisomer. The subset is drawn without replacement from the available candidate pool and capped at the pool size. The second target uses a pool of up to 1,000 candidates per molecule. These are two targets applied to existing pools, rather than independent generation runs at each target or a continuous sampling-efficiency experiment.

Both targets are defined before PoseBusters filtering. The pass mask is calculated once for the full candidate pool and projected onto the sampled subset. Rejected conformers are not replaced after filtering, so retained counts differ between methods. ChEMBL3D-PB remains the same finite stored comparison ensemble when generators receive the larger allowance. The comparisons therefore do not match retained ensemble sizes or computational cost. The in-repository classical generation procedure uses base seed 42 by default; external-pool materialization uses base seed 1729, with deterministic selections keyed by molecule and method. Executed seed overrides should be retained with the run metadata.

## Molecular identity and structural validity

For the CASF comparison, we retain conformers that pass all configured PoseBusters checks for molecular identity and intramolecular geometry. The checks cover file loading, RDKit sanitization, InChI conversion, connectivity, radicals, molecular formula, bond connectivity, tetrahedral chirality, and double-bond stereochemistry. Distance-geometry checks use configured tolerances of 0.25 for bond lengths and angles and 0.3 for clashes, excluding hydrogens. Aromatic-ring and specified trigonal double-bond planarity use a 0.25 Å flatness threshold. The energy-ratio check uses threshold 100 and a reference ensemble of 50 conformations. No protein-contact checks are included. [PoseBusters](https://doi.org/10.1039/D3SC04185A)

Generated CASF ensembles are checked against the matched molecular reference, normally the ChEMBL3D topology, whereas the stored ChEMBL3D comparison is checked against the CASF ligand. Passing establishes compliance with the specified checks; it does not guarantee compatibility with a receptor or uniformly favorable energies under the separate energy analysis. We report retained ensemble sizes alongside diversity and recovery. In the 23-molecule experiment, PoseBusters validity is assessed separately from the supplied-pool coverage calculations, as described below.

## Geometric comparison and diversity

We compare internal conformations using aligned heavy-atom RMSD. Hydrogens are removed and RDKit's `GetBestRMS` is used to account for equivalent atom mappings. If that calculation raises an exception, the CASF routine attempts `AlignMol` with atoms matched by their indices. An unresolved nonfinite reference RMSD causes that ensemble's analysis to fail. The frequency of index-based fallback is not recorded in the current exports and remains an audit item; not every score can therefore be assumed to use a verified symmetry-aware mapping.

Conformers are grouped using greedy clustering in their stored order. Each conformer is assigned to the nearest existing representative when the aligned RMSD is strictly below the clustering radius; otherwise, it starts a new cluster. Representatives are retained rather than updated. If all distances to existing representatives are nonfinite, the implementation creates a new cluster. We use a radius of 1.0 Å for the main comparison and average cluster counts over entries with defined measurements. The count depends on the size and ordering of the retained pool and describes geometric variation rather than occupied energy basins.

## Energy analysis

The planned comparison uses MMFF94s to characterize energetic plausibility and dispersion in the retained ensembles. The implemented energy helper evaluates the current geometry without an additional minimization, except that it reuses a saved post-minimization energy when a conformer is marked as minimized and that value is available. Structures without usable MMFF parameters or a finite energy return an undefined value. The current summary records the number of finite energies, minimum, maximum, median, and population standard deviation within each molecule.

For the corrected analysis, comparisons between generators should first be made within each molecule, with robust summaries across molecules and explicit reporting of unavailable energies. The intended summaries include the median of per-molecule mean energies and the median of per-molecule energy standard deviations. The current summary routine does not export the per-molecule mean; it must be calculated from the per-conformer energies or added to the corrected analysis. Saved energies also require provenance checks before they are treated as a common MMFF94s calculation. Energy spread is interpreted alongside geometry and reference recovery, not as a direct count of conformational basins or a binding free energy.

[Author note: this analysis remains pending on corrected identities and final pools. Verify the hydrogen treatment and force-field variant for saved energies, implement the required mean-energy summary, and replace this planned wording with the executed procedure. The historical energy sidecar does not supply the missing corrected result.]

## Recovery of individual bound conformations

For each CASF entry, we calculate the minimum aligned heavy-atom RMSD between the retained ensemble and its experimental reference. An entry is recovered when at least one conformer has RMSD ≤ 0.75 Å. Hit@0.75 is the fraction of all mapped entries satisfying this condition. Empty, missing, or unmeasurable outputs remain in the denominator as failures. The mean minimum RMSD is calculated only over entries with defined measurements, with that denominator reported separately. Mean retained counts use all mapped entries, including zero-count outputs.

The fraction of retained conformers close to the reference addresses how frequently a method samples near that observation, rather than whether it finds at least one match. The shared metric helper can calculate this fraction over finite RMSDs, but it is not exported by the CASF reference-scoring routine used for the current Results tables.

[Author note: if the conformer-level fraction is retained in Results, calculate it for the final valid pools and specify the handling of failed RMSDs and empty ensembles before reporting a molecule-averaged value. It must not be substituted for ligand-level recovery.]

## Molecular size and flexibility

We stratify the core panel using descriptors calculated from the experimental ligand after removing hydrogens. Heavy-atom groups are fewer than 20, 20–29, 30–39, and at least 40 atoms. Flexibility uses RDKit's `CalcNumRotatableBonds`, grouped as 0–3, 4–6, 7–8, and at least 9. The same reference-derived group assignment is used for every generator. Recovery within each group includes all mapped entries, with missing outputs counted as failures, and is accompanied by recovered and total counts. Size and flexibility are analyzed separately; neither analysis adjusts for the other descriptor or establishes an independent cause of difficulty. Groups with fewer than five entries are explicitly marked.

## Coverage of multiple experimental conformations

For each of the 23 drug-like molecules, we compare one supplied generated ensemble with all experimental references using aligned heavy-atom RMSD. COV-R is the fraction of references with at least one generated conformer at RMSD < 0.75 Å. COV-P is the fraction of generated conformers with at least one reference below the same threshold. Unlike CASF recovery, the existing coverage implementation uses a strict inequality. MAT-R averages the distance from each reference to its nearest generated conformer; MAT-P averages the distance from each generated conformer to its nearest reference. Each metric is calculated per molecule and then averaged equally across the 23 molecules.

The external evaluator removes hydrogens, attempts `GetBestRMS` for each reference–generated pair, and records failed comparisons as undefined. Reference rows and generated columns without any finite match are excluded from the corresponding coverage denominator. The imported Qwen evaluator instead includes undefined nearest-match distances as coverage misses. Both omit undefined nearest-match distances from MAT. These differences in failure handling remain part of the archived evaluation and must be harmonized before making a definitive method ranking.

Coverage and matching distances currently use the supplied pools without a common PoseBusters filter. The PB pass percentage in Table 5 is assessed separately and pools passing conformers across molecules, whereas COV and MAT give each molecule equal weight. Experimental records may be redundant and do not exhaust accessible conformational space, so an unmatched generated conformer is not automatically physically inaccessible. The molecule-level plot compares paired Qwen and LoQI outcomes for all 23 molecules; coincident points are grouped with their multiplicity, and the imatinib and actinonin examples are retained from the earlier analysis.

## Aggregation, uncertainty, and reproducibility

The working tables and figures use the corrected findings archived on 24 September 2026, with the size-stratum export from 28 September. The same ten CASF pipelines and the stored ChEMBL3D-PB ensemble are used throughout, and the 23-molecule comparison uses the five external learned generators and the same Qwen checkpoint. We preserve the distinction between whole-panel recovery and conditional means of defined RMSD or diversity measurements. Figure lines connecting the two candidate targets do not imply intermediate measurements, and the current working figures show descriptive summaries without confidence intervals.

[Author note: if the archived paired confidence intervals are retained in the final Results, describe their calculation here. The existing analysis uses 3,000 paired resamples of ChEMBL parent compounds for CASF, preserving all complex entries for each sampled compound, and 10,000 paired resamples of molecules for the drug panel. It uses percentile 95% intervals and a NumPy generator initialized with seed 20260924. These intervals do not account for training or generation variability, checkpoint selection, or multiple comparisons. They are not currently plotted in the working figures.]

Source records, exact Qwen identity, and input hashes are archived with the table and figure exports. No new generation, scoring, or force-field calculation is implied by rendering these summaries. The final dataset construction procedure, molecular selection, generation settings, and public release version remain to be specified after the evaluation with benchmark molecules excluded from Qwen training. Its principal intended application is pharmacophore-conditioned learning; the present benchmark does not test downstream training improvement.

## Implementation sources and remaining author checks

This source map supports revision and is not intended as manuscript prose. Inspected implementation and available documentation establish the procedures above; historical runtime overrides and unresolved provenance remain identified explicitly.

| Topic | Implementation and evidence |
| --- | --- |
| Plinder-23 identity selection | [Selection provenance](plinder23_appendix/selection_provenance.md); [verified identity list and counts](plinder23_appendix/molecules.csv); [Appendix A](appendix_working.md) |
| Molecular matching and stereoisomer resolution | [Data preparation](data_preparation.md); [identity checks](../src/casf_benchmark/chembl3d/identity.py); [reference loader](../src/casf_benchmark/chembl3d/loader.py) |
| Classical generation and validity configuration | [Generation procedures](generation_methods.md); [conformer pipelines](../src/casf_benchmark/generation/conformer_sets.py) |
| Candidate subsets and projected validity mask | [Materialization](materialization.md); [normalizer](../src/casf_benchmark/generation/normalizer.py) |
| RMSD, clustering, and energy | [Metric implementations](../src/casf_benchmark/analysis/metrics.py); [CASF scoring and reference cap](../src/casf_benchmark/cli/analyze_conformer_sets.py) |
| External COV/MAT | [Multiple-reference evaluator](../src/casf_benchmark/analysis/druglike_covmat.py) |
| Qwen COV/MAT failure handling | [Imported Qwen evaluator](/mnt/weka/vtarasov/code/3DMolGen/src/molgen3D/evaluation/utils.py) |
| FlowR initialization and cleanup | [FlowR adapter](/mnt/weka/vtarasov/code/flowr_root/flowr/gen/generate_conformers_from_smiles.py) |
| Checkpoint identity | [Run registry](../src/casf_benchmark/config/generation_runs.yaml); [earlier inference audit](/mnt/weka/mbedrosian/code/casf-benchmark-paper/jcim-overleaf/METHODS_NOTES.md) |
| Aggregation and archived paired intervals | [Corrected analysis](publication_analysis_2026_09_24.py); [size analysis](publication_size_analysis_2026_09_28.py) |
| Displayed numerical records | [Table provenance](results_tables/provenance.json); [figure provenance](results_figures/provenance.json) |
