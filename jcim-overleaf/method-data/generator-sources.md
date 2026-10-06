# Generator descriptions: evidence checked 6 October 2026

This record supports Supporting Table S8. Published descriptions identify each
method's representation and training objective; local adapters identify how the
benchmark invokes it. Current upstream defaults do not establish historical
runtime overrides. No missing checkpoint hash, decoding setting, or training
manifest has been inferred from a paper.

| Method | Published source | Evidence used in the table |
| --- | --- | --- |
| RDKit | [ETKDG paper](https://doi.org/10.1021/acs.jcim.0c00025); [official RDKit documentation](https://www.rdkit.org/docs/RDKit_Book.html#conformer-generation) | Distance geometry with experimental torsion preferences. Random-coordinate starts and MMFF94s processing describe our classical pipelines, not an assertion about universal RDKit defaults. |
| ChEMBL3D / LoQI | [Original paper](https://doi.org/10.26434/chemrxiv-2025-k4h7v); [author repository](https://github.com/isayevlab/LoQI); [diffusion configuration](https://github.com/isayevlab/LoQI/blob/main/scripts/conf/loqi/loqi.yaml) | Low-energy ChEMBL3D training, stereochemical graph conditioning, coordinate diffusion and a Gaussian prior. The local adapter uses the diffusion model and exports coordinates without an added minimization. Optional AIMNet2 optimization in current upstream tools is not assumed to have been applied. |
| Torsional Diffusion | [Original paper](https://arxiv.org/abs/2206.01729); [author repository](https://github.com/gcorso/torsional-diffusion) | GEOM-DRUGS checkpoint; RDKit local geometry and diffusion of torsions. The local adapter constructs RDKit seeds, randomizes torsions, and provides optional pre/post-MMFF flags that default to off. Particle guidance is not attributed to the evaluated configuration. |
| MCF drugs-L | [ICML paper](https://proceedings.mlr.press/v235/wang24q.html); [author repository](https://github.com/apple-aiml-research/ml-mcf) | GEOM-DRUGS large checkpoint, molecular coordinate fields conditioned on graph spectral features, Gaussian noise. The local adapter performs RDKit graph preparation, samples coordinates and rescales them without added minimization. |
| NExT-Mol DMT-L | [Original paper](https://arxiv.org/abs/2502.12638); [author repository](https://github.com/acharkq/NExT-Mol#3d-conformer-prediction) | The paper distinguishes DMT from DMT augmented with MoLlama. The local CASF adapter loads `drugs_dmt_l_e2999.ckpt` with `use_llm=False`. The relevant conformer training data are GEOM-DRUGS; the separate MoLlama ZINC-15 pretraining should not be credited to this evaluated DMT-L model. |
| FLOWR.root | [Paper, version 6](https://arxiv.org/html/2510.02578v6); [author checkpoint documentation](https://github.com/jule-c/flowr_root#checkpoints) | Published ligand pretraining includes ZINC3D, PubChem3D, Enamine REAL and OMol25. The public `flowr_root_v2.2.ckpt` name denotes a joint model; our launchers name `flowr_root_v2.2_mol.ckpt`, loaded by the ligand-only adapter. The exact lineage of those local weights remains unverified. The table therefore identifies the published ligand-pretraining recipe as context, rather than asserting that the local model received all later protein–ligand training stages. |
| Qwen | Author-provided model identity and project evaluation records | Initial Qwen 1.7B FSQ, step 47,023, with SMILES-conditioned coordinate tokens. The corpus was broadly assembled; exact composition and tokenizer reference remain author items. No public paper has been substituted for the unidentified tokenizer reference. |

## Local implementation evidence

Paths below are source locations in the shared workspace, not external papers.

- Classical pipelines: `src/casf_benchmark/generation/conformer_sets.py` in the analysis repository.
- LoQI: `/mnt/weka/mbedrosian/codex_dir/loqi/codex/loqi_casf_generate.py` and `scripts/conf/loqi/loqi.yaml` in the same checkout.
- Torsional Diffusion: `/mnt/weka/mbedrosian/codex_dir/torsional_diffusion/codex/torsional_diffusion_casf_generate.py` (`generate`, pre/post-MMFF flags).
- MCF: `/mnt/weka/mbedrosian/codex_dir/mcf_drugs_l/codex/mcf_casf_generate.py` (`molecule_from_smiles`, `McfSampler.generate`).
- NExT-Mol: `/mnt/weka/mbedrosian/codex_dir/nextmol_dmt_l/codex/nextmol_dmt_l_casf_generate.py` (checkpoint selection, `args.use_llm=False`, coordinate export); its upstream `model/diffusion_pl.py` provides coordinate sampling.
- FlowR: `/mnt/weka/vtarasov/code/flowr_root/flowr/gen/generate_conformers_from_smiles.py`, with `scripts/generate_conformers_core.sl`, `generate_conformers_ref.sl`, and `generate_conformers_druglike.sl` in the same checkout. These specify the local ligand-only checkpoint, harmonic graph inpainting, and force-field/stereochemical processing.

Exact hashes and executed overrides should still be recorded for final
reproducibility. The publication table summarizes mechanisms and documented
configuration, rather than filling unresolved run metadata with guessed values.
