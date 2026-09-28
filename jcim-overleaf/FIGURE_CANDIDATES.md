# Additional figure candidates

18 optional figures derived from the HTML report’s questions and corrected archived measurements. All are included at the end of the Overleaf draft through `figure-candidates.tex`. Comment out that input to hide the entire review appendix; remove individual figure blocks to select candidates. The four main figures are unchanged.

Gallery: [additional-figures.pdf](figures/candidates/additional-figures.pdf). Every candidate also has PDF, SVG, and PNG artwork.

| ID | Candidate | Report connection |
| --- | --- | --- |
| C01 | [How much recovery improves over the stored pool](figures/candidates/c01-paired-recovery.pdf) | HTML §§3–4; core-entries.csv |
| C02 | [Additional recoveries and lost recoveries](figures/candidates/c02-gained-lost.pdf) | HTML §3; core-entries.csv |
| C03 | [Which ligands are recovered by which methods?](figures/candidates/c03-ligand-recovery.pdf) | HTML §§3,6; core-entries.csv |
| C04 | [Pairwise complementarity between generators](figures/candidates/c04-complementarity.pdf) | HTML §§3,12; core-entries.csv |
| C05 | [Best available conformers versus typical conformers](figures/candidates/c05-best-typical.pdf) | HTML §3; recovery-table-tiers.csv |
| C06 | [Per-ligand geometric agreement in selected comparisons](figures/candidates/c06-paired-rmsd.pdf) | HTML §§3,5; core-entries.csv |
| C07 | [Where the core cohort has statistical support](figures/candidates/c07-cohort-composition.pdf) | HTML §§1,6; core-entries.csv |
| C08 | [Different summaries of geometric diversity](figures/candidates/c08-diversity-profiles.pdf) | HTML §5; core-entries.csv |
| C09 | [Physical validity and bound-conformer recovery](figures/candidates/c09-validity-recovery.pdf) | HTML §§3,9; core-entries.csv |
| C10 | [Classical raw and minimized pool comparisons](figures/candidates/c10-classical-pools.pdf) | HTML §9; core-entries.csv |
| C11 | [Coverage and precision for every drug molecule](figures/candidates/c11-drug-coverage-matrix.pdf) | HTML §7; drug_molecule_metrics.csv |
| C12 | [Coverage and concentration define different preferences](figures/candidates/c12-drug-tradeoff.pdf) | HTML §§7,12; drug_molecule_metrics.csv |
| C13 | [Similar mean recall can hide opposing molecule-level outcomes](figures/candidates/c13-drug-paired-molecules.pdf) | HTML §7; drug_molecule_metrics.csv |
| C14 | [Chemical validity varies strongly across molecules and methods](figures/candidates/c14-drug-validity.pdf) | HTML §7; drug_molecule_metrics.csv |
| C15 | [Continuous matching distances complement coverage](figures/candidates/c15-drug-matching-distance.pdf) | HTML §7; drug_molecule_metrics.csv |
| C16 | [The experimental reference collection is unevenly sampled](figures/candidates/c16-drug-reference-counts.pdf) | HTML §§1,7; drug_molecule_metrics.csv |
| C17 | [How uncertain are differences from LoQI?](figures/candidates/c17-drug-paired-uncertainty.pdf) | HTML §7; drug_molecule_metrics.csv |
| C18 | [Illustrative failures differ across constraints](figures/candidates/c18-drug-case-profiles.pdf) | HTML §7; drug_molecule_metrics.csv |

## Reproduction and interpretation

`python jcim-overleaf/build_candidate_figures.py` recreates the artwork and appendix using archived data only. Source hashes, cohort selection, and bootstrap settings are recorded in `figure-data/candidates/provenance.json`. Derived recovery matrices and confidence intervals are alongside the figures.

C01–C10 describe core94 and the ten main generation pipelines, with ChEMBL3D-PB where applicable. Continuous RMSD and diversity summaries are conditional on measurable ensembles; binary recovery counts missing outputs as misses. C11–C18 describe the six selected generators on the 23-molecule drug set; coverage is not harmonized by PoseBusters filtering and failure handling. No historical energy curves, random-K results, additional recovery cutoffs, reference-cohort results, or training-recipe comparisons were introduced.

These are alternatives as well as additions: C01/C02 emphasize uncertainty versus individual gains and losses; C03/C04 show individual outcomes versus pairwise complementarity; C11/C13/C18 show complete molecule profiles versus selected comparisons and examples. Keep the versions that serve the final story.
