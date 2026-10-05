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
