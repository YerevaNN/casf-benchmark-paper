# Appendix: working draft

## A. Selection and composition of the Plinder-23 panel

Plinder-23 is an exploratory, manually curated selection from PLINDER release **2024-06/v2**, rather than an official PLINDER split. Selection began with the release's ligand-per-system table. Records were grouped by the complete InChIKey, retaining stereoisomers with different full keys as separate identities. The starting collection contained **51,280 identities**; the physicochemical requirements in Table A1 retained **25,392 (49.5%)**.

**Table A1. Physicochemical screening criteria.** All conditions had to be satisfied; numerical intervals include their endpoints.

| Property | Requirement |
| --- | --- |
| RDKit parsing | Valid molecular SMILES |
| Carbon and elements | At least one C; only C, N, O, S, P, F, Cl, Br, I, B, Si, Se |
| Molecular weight | 150–650 Da |
| Crippen cLogP | −1 to 6 |
| Hydrogen-bond donors | ≤5 |
| Hydrogen-bond acceptors | ≤10 |
| Rotatable bonds | ≤12 |
| Topological polar surface area | ≤140 Å² |
| Heavy atoms | ≥10 |
| Absolute formal charge | ≤2 |

Surviving identities were ranked by descending number of distinct PDB entries, with PLINDER-system count as the secondary ordering. Ranking by PDB entries reduced the influence of multiple systems derived from one multimeric structure. Manual review then selected recognizable drugs and medicinal-chemistry-like protein ligands, excluding common crystallization additives, metabolites, cofactors, and other compounds outside the intended panel. The review also considered carbohydrates, glycosylation residues, amino acids, nucleotides, lipids, fluorescent probes, and predominantly covalently attached biological groups. These were judgment-based categories, not executable exclusion rules; bortezomib, for example, remained despite its reversible covalent mechanism.

The panel is not a random sample or simply the 23 most frequent PLINDER ligands. A later screen based on DrugBank-derived ATC annotations did not define this set, and no ATC-based inclusion requirement or minimum of ten distinct PDB entries should be inferred. Exact reproduction uses the saved molecular identities rather than repeating the manual selection. Subsequent exclusion of these molecules and related training molecules was a separate operation from panel construction.

**Table A2. Panel characteristics.** Each identity receives equal weight in the median.

| Characteristic | Minimum | Median | Maximum |
| --- | --- | --- | --- |
| Distinct PDB entries per identity | 9 | 23 | 91 |
| PLINDER systems per identity | 26 | 59 | 716 |
| Evaluated reference conformers per identity | 23 | 60 | 716 |
| Molecular weight (Da) | 162.2 | 357.8 | 505.6 |
| Heavy atoms | 11 | 24 | 37 |
| Rotatable bonds | 1 | 5 | 11 |
| Crippen cLogP | -0.86 | 1.90 | 5.70 |
| TPSA (Å²) | 16.1 | 83.6 | 131.2 |

There are **23 identities**, **2,433 system memberships** summed across the per-identity PLINDER counts, and **2,450 reference conformers** in the evaluated collection. Systems and evaluated conformers are different counting units and need not agree for an individual molecule. Distinct PDB counts are reported per identity; their sum is not a count of unique PDB entries across the panel. None of these counts represents a measured number of distinct protein targets or conformational states.

**Table A3. The 23 selected molecules in saved PDB-frequency order.** CCD denotes the Chemical Component Dictionary code; reference counts are those used in the archived evaluation.

| Rank | Molecule | CCD | Distinct PDBs | PLINDER systems | Evaluated references | MW (Da) | Heavy atoms | Rotatable bonds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | Staurosporine | STO | 91 | 179 | 177 | 466.5 | 35 | 2 |
| 2 | Actinonin | BB2 | 38 | 96 | 99 | 385.5 | 27 | 11 |
| 3 | Acetazolamide | AZM | 37 | 92 | 92 | 222.3 | 13 | 2 |
| 4 | Cyclothiazide | CYZ | 36 | 146 | 146 | 389.9 | 24 | 2 |
| 5 | Triclosan | TCL | 34 | 148 | 151 | 289.5 | 17 | 2 |
| 6 | Bortezomib | BO2 | 31 | 189 | 189 | 384.2 | 28 | 9 |
| 7 | Trimethoprim | TOP | 29 | 59 | 60 | 290.3 | 21 | 5 |
| 8 | Chloramphenicol | CLM | 28 | 160 | 157 | 323.1 | 20 | 6 |
| 9 | Imatinib | STI | 28 | 57 | 57 | 493.6 | 37 | 7 |
| 10 | Amprenavir | 478 | 25 | 30 | 30 | 505.6 | 35 | 11 |
| 11 | S-thalidomide | EF2 | 24 | 73 | 73 | 258.2 | 19 | 1 |
| 12 | Indomethacin | IMN | 23 | 51 | 54 | 357.8 | 25 | 4 |
| 13 | Nicotine | NCT | 20 | 76 | 76 | 162.2 | 12 | 1 |
| 14 | IBMX | IBM | 20 | 57 | 57 | 222.2 | 16 | 2 |
| 15 | Diclofenac | DIF | 20 | 50 | 55 | 296.2 | 19 | 4 |
| 16 | Colchicine | LOC | 20 | 30 | 30 | 399.4 | 29 | 5 |
| 17 | Dexamethasone | DEX | 19 | 32 | 32 | 392.5 | 28 | 2 |
| 18 | Fluconazole | TPF | 17 | 32 | 32 | 306.3 | 22 | 5 |
| 19 | Fosmidomycin | FOM | 16 | 30 | 28 | 183.1 | 11 | 5 |
| 20 | Bromosporine | BMF | 16 | 26 | 23 | 404.5 | 28 | 5 |
| 21 | 4-Hydroxytamoxifen | OHT | 15 | 74 | 81 | 387.5 | 29 | 8 |
| 22 | Trichostatin A | TSN | 15 | 30 | 35 | 302.4 | 22 | 6 |
| 23 | Pleconaril | W11 | 9 | 716 | 716 | 381.4 | 27 | 6 |

Molecular weights and rotatable-bond counts are taken from the saved shortlist. Heavy-atom counts are calculated from its SMILES. Full InChIKeys, SMILES, selection descriptors, and evaluated reference counts are provided in the companion [identity table](plinder23_appendix/molecules.csv). All 23 full InChIKeys match the archived evaluation identities. Counts in Table A3 refer to the saved PLINDER release and evaluated collection, not to the current contents of the PDB.

[Author note: selection provenance is now established. The extraction and preprocessing of the evaluated reference coordinates still require documentation, including why evaluated conformer counts differ from the saved system counts, protonation and tautomer handling, alternate locations, and duplicate or nearly identical structures. The present metrics count experimental reference records without a defined clustering of distinct bound states.]

## Appendix provenance notes

The source of the selection description is the author-supplied Plinder-23 selection provenance. The molecular identities and descriptors were read from the saved shortlist, and all per-identity PDB, system, and instance counts were checked against the release's ligand-per-system table. Reference counts were matched to the archived Qwen evaluation by full InChIKey. The same counts are used for the other methods in the working multi-reference comparison. No manual selection was repeated and no conformers were generated or rescored.

Source paths and input hashes are recorded in [appendix provenance](plinder23_appendix/provenance.json). The accompanying [source selection note](plinder23_appendix/selection_provenance.md) retains the supplied selection account without the surrounding conversation.


## B. Evaluation panels and generation methods

<!-- results-figure-1:start -->
![Figure S1](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-1-datasets.png)

**Figure S1. Experimental evaluation panels.** CASF-2016 is matched to ChEMBL3D by molecular identity and stereochemistry, followed by eligibility filtering, to obtain 94 core entries. The stored ChEMBL3D ensemble is the computed comparison resource. The separate 23-molecule PLINDER panel assesses several references per molecule. The larger 1,236-entry collection remains in Supporting Information. Arrows describe evaluation design, not a quantified exclusion funnel. For PLINDER 2024-06/v2, physicochemical screening reduced 51,280 identities to 25,392; ranking by distinct PDB entries and manual review yielded the 23-molecule panel (Methods; Appendix A). The supplied-pool comparison does not yet use the common CASF validity filter; reference-coordinate preprocessing remains to be documented.

[PDF](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-1-datasets.pdf) · [Editable SVG](/mnt/weka/mbedrosian/code/casf-benchmark/docs/results_figures/figure-1-datasets.svg)
<!-- results-figure-1:end -->

<!-- results-table-1:start -->
**Table S8. Evaluated generation pipelines and stored comparison ensemble.**

| Method | Representation / sampling | Training data and model | Input / initialization | Processing before common CASF filtering |
| --- | --- | --- | --- | --- |
| ChEMBL3D-PB | Stored optimized conformers | ChEMBL3D reference resource | Matched stereoisomer and its stored coordinates | PoseBusters filtering; no new generation |
| RDKit raw | ETKDGv3 distance geometry with experimental torsion preferences | No learned model | Molecular graph; random-coordinate embedding | Chirality enforcement and finite-coordinate check; no added minimization |
| RDKit minimized | ETKDGv3 distance geometry with experimental torsion preferences | No learned model | Same embedding procedure; separately sampled pool | MMFF94s minimization, up to 500 iterations |
| Torsion raw | Uniform perturbation of eligible torsions | No learned model | Stored matched geometry; displacements within ±120° | Steric clash filtering |
| Torsion minimized | Uniform perturbation of eligible torsions | No learned model | Same seed procedure; separately sampled pool | Clash filtering, MMFF94s minimization, repeat clash check |
| LoQI | Stereochemistry-aware equivariant coordinate diffusion | Low-energy ChEMBL3D structures; released diffusion model | Stereochemical molecular graph; Gaussian coordinate prior | Direct coordinate output; no added minimization in the inspected CASF adapter |
| Torsional Diffusion | Diffusion over rotatable-bond torsion angles | GEOM-DRUGS; drugs_default model | RDKit seed conformers with randomized torsions | Preserves seed local geometry; optional MMFF relaxation is disabled by default in the adapter |
| MCF drugs-L | Diffusion of coordinate fields over graph Laplacian features | GEOM-DRUGS; large MCF model | Graph spectral features and Gaussian coordinate noise; RDKit graph preparation | Coordinate rescaling; no added minimization in the inspected CASF adapter |
| NExT-Mol DMT-L | Coordinate diffusion transformer with atom and pair representations | GEOM-DRUGS; DMT-L, without MoLlama conditioning | Molecular graph; Gaussian coordinate noise | Direct coordinate output; no added minimization in the inspected CASF adapter |
| FlowR | Equivariant flow matching with the molecular graph held fixed | Published ligand pretraining: ZINC3D, PubChem3D, Enamine REAL, OMol25; local v2.2_mol model | Ligand-only generation; harmonic graph-based coordinate prior | Force-field cleanup (MMFF, UFF fallback) and stereochemistry checks |
| Qwen 1.7B FSQ | Autoregressive generation of discrete coordinate tokens | Project-assembled conformer corpus; 1.7B FSQ model, step 47,023 | Molecular SMILES; sequential coordinate-token sampling | FSQ coordinate decoding before common validity filtering |

Descriptions combine the original method papers and repositories with the inspected CASF adapters. FlowR training resources describe the published ligand-pretraining recipe; the exact training lineage of the local v2.2_mol checkpoint still requires confirmation. NExT-Mol DMT-L is evaluated without MoLlama conditioning. Optional upstream processing is not assumed to have been enabled. All CASF ensembles subsequently undergo the specified PoseBusters filtering. Exact historical checkpoint hashes and runtime overrides remain author audit items. Qwen training composition and its coordinate-tokenizer citation remain to be documented; its benchmark-excluded evaluation is pending.

Sources: [method papers, repositories, and local adapter evidence](results_tables/generator_sources.md).
<!-- results-table-1:end -->

## C. Energy summaries

**Table S9. Energy summaries at the ChEMBL-count target (kcal/mol).**

| Method | Molecules | Median energy SD | Mean energy SD | Median mean energy | Median paired mean-energy difference |
| --- | --- | --- | --- | --- | --- |
| RDKit minimized | 94 | 1.9 | 2.6 | 26.4 | -6.1 |
| Torsion minimized | 94 | 2.1 | 2.4 | 26.5 | -5.7 |
| LoQI | 94 | 1.9 | 2.8 | 34.1 | +0.2 |
| MCF drugs-L | 93 | 3.6 | 9.5 | 43.3 | +5.1 |
| NExT-Mol DMT-L | 92 | 3.0 | 4.0 | 37.9 | +3.1 |
| RDKit raw | 94 | 7.2 | 8.1 | 60.0 | +28.3 |
| Torsional Diffusion | 94 | 14.0 | 271.2 | 71.6 | +32.2 |
| FlowR | 94 | 8.8 | 12.8 | 33.5 | -0.8 |
| Torsion raw | 94 | 15.3 | 15.2 | 51.3 | +23.1 |
| Qwen 1.7B FSQ | 94 | 9.1 | 10.2 | 51.6 | +13.6 |
| ChEMBL3D-PB | 94 | 1.8 | 2.3 | 32.6 | +0.0 |

Energies use the common hydrogen preparation described in Methods. SD and mean are calculated within each molecule before aggregation; paired differences use that molecule’s ChEMBL3D mean energy. All retained conformers have finite energies and converged hydrogen relaxations; no heavy atom moved. Empty ensembles have no defined energy distribution. See the [energy protocol and records](energy_analysis/README.md).
