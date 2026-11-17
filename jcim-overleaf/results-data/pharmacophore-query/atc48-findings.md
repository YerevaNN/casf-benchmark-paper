# PLINDER-ATC48 extension panel pharmacophore-query findings

## Scope

This pilot measures same-molecule structural coverage: whether cached conformer ensembles generated from molecular identity contain rigid arrangements of direct protein-contact pharmacophore features observed in experimental complexes. It is not a virtual-screening enrichment experiment.

## Reference eligibility

- Retained references: 379
- Eligible positional queries: 311
- Eligible molecules: 44
- Query eligibility was frozen from reference-only checks before generator matching.

| Reference-only exclusion | Count |
|---|---:|
| fewer_than_four_distinct_contact_centers | 65 |
| No directly contacting protein chain | 3 |

The retained context includes nearby non-protein species for provenance, but the extractor received only directly contacting protein chains. Queries with nearby metals: 43; nearby nucleic-acid chains: 11; nearby cofactors: 96.

## Matching completeness

Pharmit supplies the pinned feature definitions and recognizes every candidate feature set. Final membership does not use its indexed triplet search: the pilot exhaustively enumerates every injective type-compatible correspondence, applies only necessary pair-distance pruning, then requires a proper rigid Kabsch fit with every feature inside the stated radius. Candidate-level membership had zero set-inclusion violations as radii increased from 0.5 to 0.75 to 1.0 Å.

## Primary endpoint (0.5 Å, budget 1,000)

| Library | Coverage | 95% molecule-bootstrap interval |
|---|---:|---:|
| Qwen/LoQI 50:50 | 0.7169 | [0.5991, 0.8289] |
| Qwen 1.7B FSQ | 0.7031 | [0.5806, 0.8147] |
| LoQI | 0.6830 | [0.5608, 0.8011] |
| RDKit ETKDG + MMFF94s | 0.5729 | [0.4553, 0.6852] |
| RDKit ETKDG raw | 0.4666 | [0.3454, 0.5882] |

## Principal paired comparison

The equal Qwen/LoQI mixture differs from Qwen alone by +0.0138 coverage points (95% interval [-0.0179, +0.0550]).

| Comparison | Difference | 95% paired-molecule interval | Expected queries gained | Expected queries lost |
|---|---:|---:|---:|---:|
| Qwen/LoQI 50:50 − Qwen 1.7B FSQ | +0.0138 | [-0.0179, +0.0550] | 10.41 | 7.02 |
| Qwen 1.7B FSQ − RDKit ETKDG raw | +0.2364 | [+0.1387, +0.3365] | 82.00 | 4.00 |
| Qwen 1.7B FSQ − RDKit ETKDG + MMFF94s | +0.1301 | [+0.0421, +0.2331] | 57.00 | 10.00 |
| LoQI − RDKit ETKDG raw | +0.2163 | [+0.1153, +0.3151] | 80.00 | 8.00 |
| LoQI − RDKit ETKDG + MMFF94s | +0.1100 | [+0.0168, +0.2073] | 52.00 | 11.00 |
| RDKit ETKDG raw − RDKit ETKDG + MMFF94s | -0.1063 | [-0.1665, -0.0504] | 11.00 | 42.00 |
| RDKit ETKDG + MMFF94s − RDKit ETKDG raw | +0.1063 | [+0.0499, +0.1641] | 42.00 | 11.00 |
| Qwen/LoQI 50:50 − RDKit ETKDG raw | +0.2502 | [+0.1562, +0.3434] | 83.11 | 1.72 |
| Qwen/LoQI 50:50 − RDKit ETKDG + MMFF94s | +0.1440 | [+0.0571, +0.2333] | 57.15 | 6.75 |

Each method uses its supplied pool without replacement; when a requested budget exceeds availability, all available conformers are used and `realized_draws` records the cap. The 1,000-conformer mixture allocates up to 500 Qwen and 500 LoQI draws. Expected gains/losses use the natural shared-Qwen coupling for mixture versus Qwen and deterministic full-pool endpoints for the single libraries.

## Prespecified sensitivities

| Radius (Å) | Library | Coverage |
|---:|---|---:|
| 0.5 | LoQI | 0.6830 |
| 0.5 | Qwen 1.7B FSQ | 0.7031 |
| 0.5 | Qwen/LoQI 50:50 | 0.7169 |
| 0.5 | RDKit ETKDG + MMFF94s | 0.5729 |
| 0.5 | RDKit ETKDG raw | 0.4666 |
| 0.75 | LoQI | 0.8155 |
| 0.75 | Qwen 1.7B FSQ | 0.8584 |
| 0.75 | Qwen/LoQI 50:50 | 0.8636 |
| 0.75 | RDKit ETKDG + MMFF94s | 0.7611 |
| 0.75 | RDKit ETKDG raw | 0.6565 |
| 1 | LoQI | 0.8916 |
| 1 | Qwen 1.7B FSQ | 0.9038 |
| 1 | Qwen/LoQI 50:50 | 0.9103 |
| 1 | RDKit ETKDG + MMFF94s | 0.8730 |
| 1 | RDKit ETKDG raw | 0.7733 |

### PoseBusters-gated sensitivity

Failed candidates remain in the sampling denominator and cannot match; no replacement or post-gate replenishment is performed.

| Library | Gated coverage |
|---|---:|
| Qwen/LoQI 50:50 | 0.7168 |
| Qwen 1.7B FSQ | 0.7031 |
| LoQI | 0.6830 |
| RDKit ETKDG + MMFF94s | 0.5729 |
| RDKit ETKDG raw | 0.4666 |

## Whole-molecule RMSD disagreement

Candidate correspondence was verified explicitly. Whole-molecule hits use RMSD ≤ 0.75 Å; pharmacophore hits use the 0.5 Å positional query.

| Library | Both | Pharmacophore only | RMSD only | Neither |
|---|---:|---:|---:|---:|
| LoQI | 47,534 | 15,401 | 35,597 | 212,468 |
| Qwen 1.7B FSQ | 41,794 | 14,518 | 29,220 | 225,112 |
| RDKit ETKDG + MMFF94s | 40,944 | 12,586 | 37,687 | 219,783 |
| RDKit ETKDG raw | 24,997 | 11,558 | 38,119 | 236,326 |

## Budget semantics and interpretation limits

All budgets count accepted cached conformers. Failed PoseBusters candidates, where that sensitivity is shown, still consume the draw and are not replenished. Generation attempts and compute are not equated across methods.

The 379 complexes are clustered within 44 selected molecules and are not independent samples. Uncertainty therefore resamples molecules. Feature definitions use the pinned local Pharmit backend and rigid matching uses the exhaustive positional correspondence check described above; no conformer was embedded, optimized, energy-ranked, or selected using the experimental references. Results measure experimental feature-arrangement coverage on this selected panel, not virtual-screening enrichment or independent model generalization.
