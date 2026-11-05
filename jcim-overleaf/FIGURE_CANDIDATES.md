# Supporting figure selection

The 18 optional candidates have been curated into Supporting Figures S3–S6.
The main manuscript no longer includes the review gallery. Supporting Figure S1
(the dataset diagram) and S2 (the approved Qwen scatter) remain.

| Kept and refreshed | Extra information beyond the main text |
| --- | --- |
| C02 → S3 | Ligands recovered only by a generator or only by ChEMBL3D; now core94 and ref1044 |
| C07 → S4 | Joint size/flexibility counts, rather than separate recovery strata; now both CASF panels |
| C11 → S5 | Coverage and precision for each of the 22 molecules and all eleven Table 4 methods |
| C14 → S6 | Validity for individual molecules, including the independent torsion baselines and four v2 Qwen models |

All other candidate figures are excluded from the compiled manuscript and SI.
The historical gallery and source records remain archived, not current results.
Repeated supporting tables were also removed; SI now has Tables S1–S6.

Current selection, exact data scope, and reproduction:
[appendix-data/README.md](appendix-data/README.md).
Run `python jcim-overleaf/appendix-data/build.py` to refresh the selected figures.
`build_candidate_figures.py` is the historical renderer and does not update the
current SI.
