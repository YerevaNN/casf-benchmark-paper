# Current ref panel: 1,044 crystals

## Primary Qwen checkpoint update, 4 November 2026

All current primary Qwen results use `qwen_1p7b_v2_coordtoken_step41054`, displayed
as **Qwen 1.7B FSQ**. It replaces step47023 in core94, ref1044 and PLINDER22;
the separate Qwen 0.6B CASF-dedup result remains. Tables 1–4, S1–S4/S8/S9,
Figures 1–6/B1 and the optional review figures use the new values. Rankings
exclude the stored ChEMBL3D baseline as before. LoQI now leads PLINDER COV-R
at 0.75 Å; Qwen retains the lowest MAT-R. Training exclusions still require
an executed manifest and must not be inferred from the run name.

Energy and sampling were not published for this checkpoint, so they were measured
from its saved PB-passing conformers under the existing paper protocols. No
coordinate optimization or validity replay was performed. The current energy
totals are 78,058/78,058 core and 734,951/734,952 ref; the one undefined ref
energy is still FlowR's. Ref has 43 empty method–entry energy pools.

See `results-data/qwen-v2/README.md` for exact source identities, reproduction
and validation. This update supersedes older checkpoint statements below;
dated records retain their original lineage. Do not restore old Qwen values
by rerunning a historical exporter.


The published cohort key is `ref_dedup`; the manuscript uses the name **ref**.
The panel contains one crystal for each of 1,044 heavy-atom isomeric SMILES
identities (1,042 ChEMBL parent IDs). Repeated crystals for a ligand do not
receive extra weight. There remain 17 parent IDs shared with core; this is
not a training-overlap exclusion or a claim of disjoint chemical parents.

The fixed seed is `casf-ref-dedup-v1-20261030`. Within each identity, selection
uses the smallest SHA-256 rank of seed, identity and crystal ID, with crystal
ID as a tie-breaker. It does not read model outputs or RMSDs. `selection.csv`
records all original candidates and which crystal was selected; `selection.json`
records the upstream algorithm and hashes. `per-entry.csv` contains only the
selected entries, the original ten paper methods at both targets, and the
stored baseline (21 × 1,044 rows).

`provenance.json` records the published SQLite hash, source cohort key, protocol
metadata and counts. `changes-from-1236.csv` compares the original and current
Tables 1–3. `counts.json` records current denominators and energy completeness.
Uncertainty is not estimated by varying this frozen selection for the paper.

## Reproduction

From the paper repository, with the benchmark checkout and published database:

```bash
python jcim-overleaf/results-data/ref-panel/export.py --database /path/to/casf_analysis_dashboard.sqlite --manifest /path/to/casf16_ref_dedup_selection_v1.csv
```

The manifest's adjacent JSON is required. This reads the published `ref_dedup`
rows and writes current ref values to the existing table, energy, sampling,
and size-stratum archives. Core and PLINDER numerical strings are preserved.
No structural calculations, force-field evaluations, or sampling simulations
are rerun. Missing RMSD is a recovery miss; continuous metrics use defined
entries; mean retained counts include zeros. Size and rotor bins use the common
ChEMBL3D-PB reference descriptors for all methods. Energy and sampling aggregates
are imported from the published common-protocol tables.

Offline rendering from the archived inputs:

```bash
python jcim-overleaf/results-data/three-dataset-tables/build_tables.py
python jcim-overleaf/results-data/ref-panel/build_supporting.py
python jcim-overleaf/results-data/core-ref-protocol/build_energy_table.py
python jcim-overleaf/results-data/three-dataset-tables/plot_casf.py
python jcim-overleaf/results-data/core-ref-protocol/build_figures.py
python jcim-overleaf/results-data/size-flexibility/build_figure.py
python jcim-overleaf/plinder22-data/plot_datasets.py
python jcim-overleaf/results-data/ref-panel/verify.py
python jcim-overleaf/plinder22-data/verify.py
```

Updated displays: Tables 1–3, S2/S3/S9, Figures 2–5/B1, and Figure S1.
Figure 5 remains CASF-only. All PLINDER figures and Table 4 retain their scope
and values. Historical ref1236 exports remain only as dated provenance; they
are not inputs to the active renderers above. The full manuscript is not
compiled locally; review it in Overleaf.
