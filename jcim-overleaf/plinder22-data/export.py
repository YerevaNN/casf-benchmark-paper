"""Archive the cleaned PLINDER-22 panel and validate it against its RMSD caches.

Run in the benchmark environment: python export.py /path/to/casf-benchmark
No generation, alignment, clustering, or energy calculation is performed.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from casf_benchmark.analysis.druglike_covmat import load_rmsd_matrices
from casf_benchmark.analysis.plinder_panel import DEFAULT_REFERENCE_DIR, load_reference_sdf_dir

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
BENCHMARK = Path(sys.argv[1]).resolve()
METHODS = {
    'rdkit_random_raw_druglike': 'RDKit raw',
    'torsion_raw_druglike': 'Torsion raw',
    'loqi_druglike': 'LoQI',
    'torsional_diffusion_druglike': 'Torsional Diffusion',
    'mcf_drugs_l_druglike': 'MCF drugs-L',
    'nextmol_dmt_l_druglike': 'NExT-Mol DMT-L',
    'flowr_druglike': 'FlowR',
    'qwen_1p7b_fsq_bigdata_step47023': 'Qwen 1.7B FSQ',
}
references = load_reference_sdf_dir(DEFAULT_REFERENCE_DIR)
assert len(references) == 22
assert sum(len(r['confs']) for r in references.values()) == 1067
original_identities = pd.read_csv(PAPER / 'method-data/plinder23-molecules.csv').set_index('ccd')
for smiles, row in references.items():
    assert Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)) == original_identities.loc[row['ccd'], 'inchikey']

sources = {}
def record(path):
    sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path
for path in sorted(DEFAULT_REFERENCE_DIR.glob('*.sdf')):
    record(path)
for name in ['analysis/plinder_panel.py', 'analysis/druglike_covmat.py', 'cli/analyze_plinder.py']:
    record(BENCHMARK / 'src/casf_benchmark' / name)
record(DEFAULT_REFERENCE_DIR.parents[1] / 'README.md')
(PAPER / 'plinder22-data/reference-source-notes.md').write_text((DEFAULT_REFERENCE_DIR.parents[1] / 'README.md').read_text())

inventory = []
for smiles, row in references.items():
    mols = row['confs']
    inventory.append(dict(smiles=smiles, name=row['names'][0], ccd=row['ccd'],
        num_true_confs=len(mols), n_pdbs=len({m.GetProp('pdb_id') for m in mols}),
        n_systems=len({m.GetProp('system_id') for m in mols}),
        represented_raw_references=sum(int(m.GetProp('cluster_size')) for m in mols)))
pd.DataFrame(inventory).to_csv(HERE / 'reference-inventory.csv', index=False)
records = []
for smiles, row in references.items():
    for mol in row['confs']:
        records.append({'smiles': smiles, **{key: mol.GetProp(key) for key in mol.GetPropNames()}})
pd.DataFrame(records).to_csv(HERE / 'reference-records.csv', index=False)
assert len({r['pdb_id'] for r in records}) == 552
assert len({r['system_id'] for r in records}) == 1045
assert sum(int(r['cluster_size']) for r in records) == 2258


frames, summaries, curve_rows, entries, minima = [], [], [], [], {}
thresholds = np.linspace(0, 1.5, 301)
for mi, (method, display) in enumerate(METHODS.items()):
    folder = BENCHMARK / 'data/results/plinder_panel' / method
    frame = pd.read_csv(record(folder / 'per_molecule.csv'))
    summary = pd.read_csv(record(folder / 'summary.csv')).iloc[0]
    cached_curves = pd.read_csv(record(folder / 'recovery.csv'))
    cache = record(BENCHMARK / 'data/results/cache/plinder_rmsd_matrices' / f'{method}.npz')
    matrices = load_rmsd_matrices(cache)
    assert len(frame) == frame.smiles.nunique() == 22
    assert set(frame.smiles) == set(matrices) == set(references)
    frame = frame.set_index('smiles').loc[list(references)].reset_index()
    recalls, precisions = [], []
    for ei, row in enumerate(frame.itertuples()):
        matrix = np.asarray(matrices[row.smiles], dtype=float)
        assert matrix.shape == (len(references[row.smiles]['confs']), row.n_generated)
        matrix = np.where(np.isfinite(matrix), matrix, np.inf)
        r, p = matrix.min(axis=1), matrix.min(axis=0)
        assert np.isfinite(r).all() and np.isfinite(p).all()
        assert (r >= 0).all() and (p >= 0).all()
        np.testing.assert_allclose([r.mean(), p.mean()], [row.mat_r, row.mat_p], rtol=0, atol=1e-12)
        recall = (r[:, None] <= thresholds).mean(axis=0)
        precision = (p[:, None] <= thresholds).mean(axis=0)
        for threshold, tag in [(0.5, '0p5'), (0.75, '0p75'), (1., '1p0')]:
            index = int(round(threshold / .005))
            np.testing.assert_allclose([recall[index], precision[index]],
                [getattr(row, 'cov_r_' + tag), getattr(row, 'cov_p_' + tag)], rtol=0, atol=1e-12)
        assert row.conformer_count == row.n_generated
        assert 0 <= row.n_pb_pass <= row.n_generated
        key = f'm{mi}_e{ei}'
        entries.append(dict(method=method, smiles=row.smiles, key=key))
        minima[key + '_reference'], minima[key + '_generated'] = r, p
        recalls.append(recall); precisions.append(precision)
    recall, precision = np.mean(recalls, axis=0), np.mean(precisions, axis=0)
    for values in [recall, precision]:
        assert np.all(np.diff(values) >= -1e-12) and np.all((values >= 0) & (values <= 1))
    for row in cached_curves.itertuples():
        index = int(round(row.threshold / .005))
        np.testing.assert_allclose([recall[index], precision[index]], [row.cov_r_mean, row.cov_p_mean], rtol=0, atol=1e-12)
    for source, dest in [('greedy_clusters_0p5','mean_clusters_0p5'),
                         ('greedy_clusters_1p0','mean_clusters_1p0'),
                         ('largest_cluster_fraction_1p0','largest_cluster_fraction_1p0'),
                         ('mat_r','mat_r_mean'), ('mat_p','mat_p_mean')]:
        np.testing.assert_allclose(frame[source].mean(), summary[dest], rtol=0, atol=1e-12)
    frames.append(frame); summaries.append(summary)
    curve_rows.extend(dict(method=method, threshold=t, recall=r, precision=p)
                      for t, r, p in zip(thresholds, recall, precision))
all_rows = pd.concat(frames, ignore_index=True)
all_rows.to_csv(HERE / 'per-molecule.csv', index=False)
pd.DataFrame(summaries).to_csv(HERE / 'summary.csv', index=False)
pd.DataFrame(curve_rows).to_csv(HERE / 'threshold-curves.csv', index=False)
np.savez_compressed(HERE / 'nearest-distances.npz', **minima)
(HERE / 'entries.json').write_text(json.dumps(entries, indent=2) + '\n')
(HERE / 'methods.json').write_text(json.dumps(METHODS, indent=2) + '\n')
# One consistent source for Tables 2/4/S4 and candidate figures.
all_rows.rename(columns={'cov_r_0p75':'cov_r_075','cov_p_0p75':'cov_p_075'}).to_csv(
    PAPER / 'results-data/druglike_selected.csv', index=False)
# Preserve every CASF cell and the approved wide table order.
for number in [1, 2]:
    path = PAPER / f'results-data/three-dataset-tables/table{number}.csv'
    with path.open() as f:
        reader = csv.DictReader(f); fields = reader.fieldnames; rows = list(reader)
    by_name = {METHODS[s['label']]: s for s in summaries}
    for row in rows:
        if row['cohort'] != 'plinder':
            continue
        row['n_total'] = 22
        s = by_name.get(row['method'])
        if s is None:
            continue
        if number == 1:
            row.update(n_defined=22, clusters_0p5=s.mean_clusters_0p5,
                clusters_1p0=s.mean_clusters_1p0,
                largest_percent=100*s.largest_cluster_fraction_1p0,
                status='PLINDER-22 full supplied pool; cleaned references 2026-10-20')
        else:
            row.update(n_distance=22, recovery_0p5=100*s.cov_r_0p5_mean,
                recovery_0p75=100*s.cov_r_0p75_mean, distance=s.mat_r_mean)
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
(HERE / 'provenance.json').write_text(json.dumps(dict(
    release='Cleaned PLINDER references, 2026-10-20', molecules=22, references=1067,
    source_sha256=sources, methods=METHODS,
    export_environment=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__, rdkit=rdBase.rdkitVersion),
    clustering='Full supplied pools; no PB filter; greedy radii 0.5/1.0 angstrom',
    coverage='RMSD <= threshold; all reference and generated minima finite; equal molecule weights',
    pb='Separate geometry/identity evaluation; not applied to coverage or clustering',
    qwen='Original 1.7B FSQ step47023; no retrained checkpoint substituted',
    reference_weight='Equal deduplicated representatives within each molecule; not weighted by cluster_size',
    checks=['full InChIKeys match the original selection', '22 identical molecular identities for all eight methods', '1067 references per method',
            'cache shapes agree with full generated pools', 'all nearest distances finite',
            'per-molecule and aggregate COV/MAT match released CSVs within 1e-12',
            'published recovery curves agree within 1e-12', 'all curves bounded and monotonic',
            'cluster counts and occupancy reproduce released summaries', 'CASF CSV cells preserved']
), indent=2) + '\n')
print('Validated 8 methods × 22 molecules against 1,067 cleaned references; archived unified PLINDER sources.')
