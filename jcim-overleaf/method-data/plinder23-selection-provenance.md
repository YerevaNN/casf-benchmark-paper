Plinder-23 Selection Provenance
Status and intended use
Plinder-23 is an internal exploratory set of 23 ligand identities selected from the PLINDER 2024-06/v2 release. It is not an official PLINDER split or a canonical benchmark. The set was later used both for bioactive-conformer evaluation and as a query set when removing related molecules from FSQ and coordinate-bin training corpora.

The deterministic physicochemical screen and structure-frequency ranking can be reproduced. The final choice of 23 molecules cannot be regenerated from those rules alone because it included a manual semantic-curation step. To reproduce the existing set exactly, reuse the identities in druglike_curated_shortlist.tsv.

Source data and molecular identity
The selection began from PLINDER’s 2024-06/v2/fingerprints/ligands_per_system.parquet. The working copy and derived tables are under:

/mnt/weka/hrant/plinder_frequency/2024-06-v2
Rows were grouped by the complete inchikeys value. Consequently, stereoisomers with different full InChIKeys remained separate. The following coverage counts were calculated for each identity:

count(DISTINCT system_id) AS n_systems
count(*)                  AS n_instances
count(DISTINCT pdb_id)    AS n_pdbs
The starting table contained 51,280 standardized ligand identities.

Deterministic physicochemical screen
SMILES were parsed with RDKit. An identity survived only when every condition below held:

Property	Requirement
RDKit parsing	MolFromSmiles(smiles) != None
Carbon	At least one carbon atom
Elements	Subset of C, N, O, S, P, F, Cl, Br, I, B, Si, Se
Molecular weight	150–650 Da
RDKit Crippen cLogP	-1 to 6
Hydrogen-bond donors	At most 5
Hydrogen-bond acceptors	At most 10
Rotatable bonds	At most 12
TPSA	At most 140 A²
Heavy atoms	At least 10
Absolute formal charge	At most 2
The implementation-equivalent predicate was:

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski

mol = Chem.MolFromSmiles(smiles)
allowed = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "Se"}

if mol is None:
    keep = False
else:
    elements = {atom.GetSymbol() for atom in mol.GetAtoms()}
    keep = (
        elements <= allowed
        and "C" in elements
        and 150 <= Descriptors.MolWt(mol) <= 650
        and -1 <= Crippen.MolLogP(mol) <= 6
        and Lipinski.NumHDonors(mol) <= 5
        and Lipinski.NumHAcceptors(mol) <= 10
        and Lipinski.NumRotatableBonds(mol) <= 12
        and Descriptors.TPSA(mol) <= 140
        and Lipinski.HeavyAtomCount(mol) >= 10
        and abs(Chem.GetFormalCharge(mol)) <= 2
    )
This stage retained 25,392 identities. Its output is:

/mnt/weka/hrant/plinder_frequency/2024-06-v2/druglike_physchem_screen.tsv
Ranking and manual selection
Survivors were ranked primarily by descending n_pdbs, with PLINDER-system count used as the effective secondary ordering. Distinct PDB entries were preferred to raw system count because one multimeric PDB structure can generate many PLINDER systems.

The final 23 were then selected manually as recognizable drugs or medicinal-chemistry-like protein ligands. The manual review excluded obvious ions and inorganic compounds; crystallization solvents, buffers, and detergents; carbohydrates and glycosylation residues; amino acids and common endogenous metabolites; nucleotides and common cofactors; predominantly covalently attached biological groups; lipids and membrane additives; and fluorescent probes such as ANS.

These categories describe the judgment used at the time, not executable exclusion rules. They were not applied with perfect categorical consistency; for example, bortezomib remained in the shortlist despite its reversible covalent mechanism. This is another reason to treat the set as exploratory.

The resulting set, in its saved PDB-frequency order, is:

Rank	Molecule	CCD	Distinct PDBs	PLINDER systems
1	Staurosporine	STO	91	179
2	Actinonin	BB2	38	96
3	Acetazolamide	AZM	37	92
4	Cyclothiazide	CYZ	36	146
5	Triclosan	TCL	34	148
6	Bortezomib	BO2	31	189
7	Trimethoprim	TOP	29	59
8	Chloramphenicol	CLM	28	160
9	Imatinib	STI	28	57
10	Amprenavir	478	25	30
11	S-thalidomide	EF2	24	73
12	Indomethacin	IMN	23	51
13	Nicotine	NCT	20	76
14	IBMX	IBM	20	57
15	Diclofenac	DIF	20	50
16	Colchicine	LOC	20	30
17	Dexamethasone	DEX	19	32
18	Fluconazole	TPF	17	32
19	Fosmidomycin	FOM	16	30
20	Bromosporine	BMF	16	26
21	4-Hydroxytamoxifen	OHT	15	74
22	Trichostatin A	TSN	15	30
23	Pleconaril	W11	9	716
The authoritative identity list, including full InChIKeys, SMILES, rotatable-bond counts, and saved coverage fields, is:

/mnt/weka/hrant/plinder_frequency/2024-06-v2/druglike_curated_shortlist.tsv
Later ATC screen was not the source of Plinder-23
A later, more objective screen resolved frequent candidates against the RCSB Chemical Component Dictionary and used DrugBank-derived ATC annotations. It found 63 annotated medicines among 270 resolved identities with at least 10 PDB entries. That analysis produced:

/mnt/weka/hrant/plinder_frequency/2024-06-v2/druglike_atc_candidates.tsv
It did not retroactively define the original Plinder-23. In particular, valid medicines found by the ATC screen can be absent from the manual shortlist, while the shortlist includes exploratory bioactive compounds without an ATC-based inclusion requirement.

How the set entered training-data filtering
The filtering audit read druglike_curated_shortlist.tsv, used each row’s smiles as the query structure, and labeled it plinder23. The resulting 23 canonical isomeric SMILES appear in:

/mnt/weka/hrant/dodock_similarity_audit_20260914_casf285_plinder23_only/queries.csv
Those 23 identities and their strict Morgan/Tanimoto > 0.4 neighborhoods were removed from the filtered training corpora. This neighborhood filtering is separate from, and subsequent to, the selection process documented here.

Reproducibility guidance
To reproduce the existing experiments, use the saved 23 identities exactly; do not rerun manual curation and assume it will yield the same list.
To build a publication-grade successor benchmark, replace manual semantic curation with explicit external annotations and exclusion rules, including a declared drug database, minimum distinct PDB/protein coverage, treatment of covalent complexes, and fixed cofactor/metabolite/additive exclusion lists.
Do not describe Plinder-23 as a random sample, the 23 most frequent PLINDER ligands, an official PLINDER split, or a fully algorithmic drug set.
