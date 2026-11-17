# pharmacophore_plinder22 pharmacophore-query findings

## Scope

This pilot measures same-molecule structural coverage: whether cached conformer ensembles generated from molecular identity contain rigid arrangements of direct protein-contact pharmacophore features observed in experimental complexes. It is not a virtual-screening enrichment experiment.

## Reference eligibility

- Retained references: 1,067
- Eligible positional queries: 932
- Eligible molecules: 22
- Query eligibility was frozen from reference-only checks before generator matching.

| Reference-only exclusion | Count |
|---|---:|
| fewer_than_four_distinct_contact_centers | 133 |
| covalent_complex_excluded | 2 |

The retained context includes nearby non-protein species for provenance, but the extractor received only directly contacting protein chains. Queries with nearby metals: 209; nearby nucleic-acid chains: 3; nearby cofactors: 414.

## Qwen completion

The immutable original pool contained 21,947 accepted conformers. The unbiased top-up made 59 attempts, accepted 53 additions, and rejected 6; every molecule then had exactly 1,000 accepted cached conformers. Original artifact hashes and original-prefix coordinate hashes were unchanged.

## Matching completeness

Pharmit supplies the pinned feature definitions and recognizes every candidate feature set. Final membership does not use its indexed triplet search: the pilot exhaustively enumerates every injective type-compatible correspondence, applies only necessary pair-distance pruning, then requires a proper rigid Kabsch fit with every feature inside the stated radius. Candidate-level membership had zero set-inclusion violations as radii increased from 0.5 to 0.75 to 1.0 Å.

## Primary endpoint (0.5 Å, budget 1,000)

| Library | Coverage | 95% molecule-bootstrap interval |
|---|---:|---:|
| Qwen/LoQI 50:50 | 0.7812 | [0.6529, 0.9001] |
| LoQI | 0.7363 | [0.5886, 0.8650] |
| Qwen 1.7B FSQ | 0.7122 | [0.5674, 0.8504] |
| Torsional Diffusion | 0.7058 | [0.5524, 0.8496] |
| MCF | 0.6993 | [0.5540, 0.8368] |
| FlowR | 0.6908 | [0.5492, 0.8266] |
| NExT-Mol | 0.6713 | [0.5147, 0.8164] |
| RDKit minimized | 0.5293 | [0.3841, 0.6706] |
| RDKit raw | 0.5211 | [0.3577, 0.6840] |

## Principal paired comparison

The equal Qwen/LoQI mixture differs from Qwen alone by +0.0690 coverage points (95% interval [+0.0026, +0.1544]).

| Comparison | Difference | 95% paired-molecule interval | Expected queries gained | Expected queries lost |
|---|---:|---:|---:|---:|
| Qwen/LoQI 50:50 − Qwen 1.7B FSQ | +0.0690 | [+0.0026, +0.1544] | 56.37 | 22.13 |
| Qwen 1.7B FSQ − RDKit raw | +0.1911 | [+0.0742, +0.3196] | 233.00 | 31.00 |
| Qwen 1.7B FSQ − RDKit minimized | +0.1829 | [+0.0642, +0.3178] | 249.00 | 31.00 |
| LoQI − RDKit raw | +0.2152 | [+0.1003, +0.3395] | 208.00 | 20.00 |
| LoQI − RDKit minimized | +0.2070 | [+0.1161, +0.3058] | 214.00 | 10.00 |
| FlowR − RDKit raw | +0.1697 | [+0.0828, +0.2679] | 178.00 | 26.00 |
| FlowR − RDKit minimized | +0.1615 | [+0.0760, +0.2636] | 187.00 | 19.00 |
| MCF − RDKit raw | +0.1782 | [+0.0803, +0.2903] | 181.00 | 25.00 |
| MCF − RDKit minimized | +0.1700 | [+0.0854, +0.2684] | 188.00 | 16.00 |
| NExT-Mol − RDKit raw | +0.1502 | [+0.0254, +0.2782] | 182.00 | 43.00 |
| NExT-Mol − RDKit minimized | +0.1421 | [+0.0382, +0.2482] | 186.00 | 31.00 |
| Torsional Diffusion − RDKit raw | +0.1847 | [+0.0934, +0.2867] | 191.00 | 12.00 |
| Torsional Diffusion − RDKit minimized | +0.1765 | [+0.0779, +0.2873] | 214.00 | 19.00 |
| RDKit raw − RDKit minimized | -0.0082 | [-0.1090, +0.0905] | 100.00 | 84.00 |
| RDKit minimized − RDKit raw | +0.0082 | [-0.0963, +0.1138] | 84.00 | 100.00 |
| Qwen/LoQI 50:50 − RDKit raw | +0.2601 | [+0.1455, +0.3872] | 247.43 | 11.19 |
| Qwen/LoQI 50:50 − RDKit minimized | +0.2519 | [+0.1527, +0.3646] | 260.11 | 7.86 |

Each method uses its supplied pool without replacement; when a requested budget exceeds availability, all available conformers are used and `realized_draws` records the cap. The 1,000-conformer mixture allocates up to 500 Qwen and 500 LoQI draws. Expected gains/losses use the natural shared-Qwen coupling for mixture versus Qwen and deterministic full-pool endpoints for the single libraries.

## Prespecified sensitivities

| Radius (Å) | Library | Coverage |
|---:|---|---:|
| 0.5 | FlowR | 0.6908 |
| 0.5 | LoQI | 0.7363 |
| 0.5 | MCF | 0.6993 |
| 0.5 | NExT-Mol | 0.6713 |
| 0.5 | Qwen 1.7B FSQ | 0.7122 |
| 0.5 | Qwen/LoQI 50:50 | 0.7812 |
| 0.5 | RDKit minimized | 0.5293 |
| 0.5 | RDKit raw | 0.5211 |
| 0.5 | Torsional Diffusion | 0.7058 |
| 0.75 | FlowR | 0.8259 |
| 0.75 | LoQI | 0.8425 |
| 0.75 | MCF | 0.8317 |
| 0.75 | NExT-Mol | 0.8198 |
| 0.75 | Qwen 1.7B FSQ | 0.8619 |
| 0.75 | Qwen/LoQI 50:50 | 0.8992 |
| 0.75 | RDKit minimized | 0.7061 |
| 0.75 | RDKit raw | 0.6870 |
| 0.75 | Torsional Diffusion | 0.7945 |
| 1 | FlowR | 0.8951 |
| 1 | LoQI | 0.9101 |
| 1 | MCF | 0.9028 |
| 1 | NExT-Mol | 0.8697 |
| 1 | Qwen 1.7B FSQ | 0.9538 |
| 1 | Qwen/LoQI 50:50 | 0.9587 |
| 1 | RDKit minimized | 0.8156 |
| 1 | RDKit raw | 0.7616 |
| 1 | Torsional Diffusion | 0.8376 |

Original Qwen-pool coverage was 0.7122; completed-Qwen coverage was 0.7122 (marginal difference +0.0000).

### PoseBusters-gated sensitivity

Failed candidates remain in the sampling denominator and cannot match; no replacement or post-gate replenishment is performed.

| Library | Gated coverage |
|---|---:|
| Qwen/LoQI 50:50 | 0.7762 |
| LoQI | 0.7298 |
| Qwen 1.7B FSQ | 0.7102 |
| Torsional Diffusion | 0.7053 |
| MCF | 0.6764 |
| NExT-Mol | 0.6433 |
| FlowR | 0.6409 |
| RDKit minimized | 0.5293 |
| RDKit raw | 0.5211 |

| Library | Archived passes | Comparable replay | Difference | Replacement passes |
|---|---:|---:|---:|---:|
| Qwen 1.7B FSQ | 20,296 | 20,297 | +1 | 50 |
| LoQI | 19,948 | 19,948 | +0 | 0 |
| RDKit raw | 22,000 | 22,000 | +0 | 0 |
| RDKit minimized | 21,999 | 21,999 | +0 | 0 |

## Whole-molecule RMSD disagreement

Candidate correspondence was verified explicitly. Whole-molecule hits use RMSD ≤ 0.75 Å; pharmacophore hits use the 0.5 Å positional query.

| Library | Both | Pharmacophore only | RMSD only | Neither |
|---|---:|---:|---:|---:|
| FlowR | 52,174 | 40,749 | 55,846 | 756,748 |
| LoQI | 160,314 | 46,645 | 109,614 | 615,427 |
| MCF | 46,442 | 59,688 | 71,657 | 754,213 |
| NExT-Mol | 48,124 | 60,202 | 57,626 | 766,048 |
| Qwen 1.7B FSQ | 45,358 | 61,227 | 69,201 | 756,214 |
| RDKit minimized | 122,894 | 58,416 | 105,985 | 644,705 |
| RDKit raw | 89,454 | 56,371 | 76,521 | 709,654 |
| Torsional Diffusion | 116,294 | 48,097 | 74,396 | 693,213 |

## Budget semantics and interpretation limits

All budgets count accepted cached conformers. Failed PoseBusters candidates, where that sensitivity is shown, still consume the draw and are not replenished. The extra Qwen generation effort is reported separately in `qwen_topup/`.

The 1,067 complexes are clustered within 22 selected molecules and are not independent samples. Uncertainty therefore resamples molecules. Feature definitions use the pinned local Pharmit backend and rigid matching uses the exhaustive positional correspondence check described above; no conformer was embedded, optimized, energy-ranked, or selected using the experimental references. Results measure experimental feature-arrangement coverage on this selected panel, not virtual-screening enrichment or independent model generalization.
