"""Import the two independent PLINDER torsion runs without changing other methods."""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import shutil

import numpy as np
import pandas as pd
from casf_benchmark.analysis.druglike_covmat import load_rmsd_matrices

HERE = Path(__file__).resolve().parent
PANEL = HERE.parent
PAPER = PANEL.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('benchmark', type=Path)
args = parser.parse_args()
benchmark = args.benchmark.resolve()
mapping = {
    'torsion_raw_druglike': 'torsion_independent_raw_druglike',
    'torsion_minimized_22_druglike': 'torsion_independent_minimized_druglike',
}
display = dict(zip(mapping.values(), ['Torsion raw', 'Torsion minimized']))
sources = {}


def record(path):
    sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def save_frame(path, frame, key='label', changed=None, only_plinder=False):
    """Preserve the original CSV cells for every observation outside this update."""
    changed = set(mapping) | set(mapping.values()) if changed is None else changed
    with path.open() as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        old = list(reader)
    new = list(csv.DictReader(io.StringIO(frame.to_csv(index=False))))
    assert len(old) == len(new)
    rows = []
    for a, b in zip(old, new):
        replace = a[key] in changed and (not only_plinder or a['cohort'] == 'plinder')
        if not replace:
            assert a[key] == b[key]
        rows.append(b if replace else a)
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def replace_rows(path, new, key='label'):
    old = pd.read_csv(path)
    labels = set(mapping) | set(mapping.values())
    before = old[~old[key].isin(labels)].reset_index(drop=True)
    chunks = []
    for method, group in old.groupby(key, sort=False):
        target = mapping.get(method, method)
        if target in new:
            replacement = new[target].copy().set_index('smiles')
            replacement = replacement.loc[group.smiles].reset_index()
            chunks.append(replacement.reindex(columns=old.columns))
        else:
            chunks.append(group)
    result = pd.concat(chunks, ignore_index=True)
    pd.testing.assert_frame_equal(before, result[~result[key].isin(labels)].reset_index(drop=True))
    save_frame(path, result, key)
    return result


inventory = pd.read_csv(PANEL / 'reference-inventory.csv').set_index('smiles')
assert len(inventory) == 22 and inventory.num_true_confs.sum() == 1067
entries = json.loads((PANEL / 'entries.json').read_text())
with np.load(PANEL / 'nearest-distances.npz', allow_pickle=False) as archive:
    minima = {key: archive[key].copy() for key in archive.files}
frames, summaries, curves = {}, {}, {}
thresholds = np.linspace(0, 1.5, 301)
for old, method in mapping.items():
    upstream = benchmark / 'data/results/plinder_panel' / method
    generation = benchmark / 'data/results/generation/plinder_independent_torsion'
    for name in ['per_molecule.csv', 'summary.csv', 'recovery.csv']:
        shutil.copy2(record(upstream / name), HERE / f'{method}-{name}')
    manifest = record(generation / method / 'generation_summary.json')
    shutil.copy2(manifest, HERE / f'{method}-generation.json')
    info = json.loads(manifest.read_text())
    assert info['experimental_coordinates_used'] is False
    assert info['seed_source'] == 'rdkit_etkdg_v3'
    assert info['n_etkdg_starts'] == 5 and info['conformers_per_start'] == 200
    assert info['n_success'] == 22 and info['total_conformers'] == 22000
    frame = pd.read_csv(upstream / 'per_molecule.csv')
    assert len(frame) == frame.smiles.nunique() == 22
    assert set(frame.smiles) == set(inventory.index)
    assert frame.n_generated.eq(1000).all() and frame.conformer_count.eq(1000).all()
    np.testing.assert_array_equal(frame.num_true_confs, inventory.loc[frame.smiles, 'num_true_confs'])
    matrices = load_rmsd_matrices(record(benchmark / 'data/results/cache/plinder_rmsd_matrices' / f'{method}.npz'))
    assert set(matrices) == set(frame.smiles)
    recall, precision = [], []
    for row in frame.itertuples():
        matrix = np.asarray(matrices[row.smiles], dtype=float)
        assert matrix.shape == (row.num_true_confs, row.n_generated)
        r, g = matrix.min(axis=1), matrix.min(axis=0)
        assert np.isfinite(r).all() and np.isfinite(g).all()
        np.testing.assert_allclose([r.mean(), g.mean()], [row.mat_r, row.mat_p], atol=1e-12, rtol=0)
        for cutoff, tag in [(.5, '0p5'), (.75, '0p75'), (1., '1p0')]:
            np.testing.assert_allclose([(r <= cutoff).mean(), (g <= cutoff).mean()],
                                      [getattr(row, 'cov_r_' + tag), getattr(row, 'cov_p_' + tag)], atol=1e-12, rtol=0)
        entry = next(e for e in entries if e['method'] in [old, method] and e['smiles'] == row.smiles)
        entry['method'] = method
        minima[entry['key'] + '_reference'] = r
        minima[entry['key'] + '_generated'] = g
        recall.append((r[:, None] <= thresholds).mean(axis=0))
        precision.append((g[:, None] <= thresholds).mean(axis=0))
    recall, precision = np.mean(recall, axis=0), np.mean(precision, axis=0)
    published = pd.read_csv(upstream / 'recovery.csv')
    for row in published.itertuples():
        index = round(row.threshold / .005)
        np.testing.assert_allclose([recall[index], precision[index]], [row.cov_r_mean, row.cov_p_mean], atol=1e-12, rtol=0)
    frames[method] = frame
    summaries[method] = pd.read_csv(upstream / 'summary.csv').iloc[0]
    for source, target in [('mat_r', 'mat_r_mean'), ('mat_p', 'mat_p_mean'),
                           ('greedy_clusters_0p5', 'mean_clusters_0p5'),
                           ('greedy_clusters_1p0', 'mean_clusters_1p0'),
                           ('largest_cluster_fraction_1p0', 'largest_cluster_fraction_1p0')]:
        np.testing.assert_allclose(frame[source].mean(), summaries[method][target], atol=1e-12, rtol=0)
    curves[method] = pd.DataFrame(dict(method=method, threshold=thresholds, recall=recall, precision=precision))

replace_rows(PANEL / 'per-molecule.csv', frames)
replace_rows(PAPER / 'results-data/druglike_selected.csv', {
    m: f.rename(columns={'cov_r_0p75': 'cov_r_075', 'cov_p_0p75': 'cov_p_075'}) for m, f in frames.items()})
replace_rows(PANEL / 'molecule-tracking/molecule_matching_values.csv', frames)
summary = pd.read_csv(PANEL / 'summary.csv')
for old, method in mapping.items():
    index = summary.index[summary.label.isin([old, method])]
    assert len(index) == 1
    summary.loc[index[0]] = summaries[method].reindex(summary.columns)
save_frame(PANEL / 'summary.csv', summary)
all_curves = pd.read_csv(PANEL / 'threshold-curves.csv')
all_curves = pd.concat([curves.get(mapping.get(m, m), g) for m, g in all_curves.groupby('method', sort=False)], ignore_index=True)
save_frame(PANEL / 'threshold-curves.csv', all_curves, key='method')
np.savez_compressed(PANEL / 'nearest-distances.npz', **minima)
(PANEL / 'entries.json').write_text(json.dumps(entries, indent=2) + '\n')
methods = json.loads((PANEL / 'methods.json').read_text())
methods = {mapping.get(m, m): label for m, label in methods.items()}
(PANEL / 'methods.json').write_text(json.dumps(methods, indent=2) + '\n')
for number in [1, 2]:
    path = PAPER / f'results-data/three-dataset-tables/table{number}.csv'
    table = pd.read_csv(path)
    for method, label in display.items():
        s = summaries[method]
        update = dict(n_defined=22, clusters_0p5=s.mean_clusters_0p5, clusters_1p0=s.mean_clusters_1p0,
                      largest_percent=100*s.largest_cluster_fraction_1p0,
                      status='PLINDER-22 independent ETKDG torsion; full supplied pools') if number == 1 else dict(
                          n_distance=22, recovery_0p5=100*s.cov_r_0p5_mean, recovery_0p75=100*s.cov_r_0p75_mean, distance=s.mat_r_mean)
        mask = table.cohort.eq('plinder') & table.method.eq(label)
        assert mask.sum() == 1
        for column, value in update.items(): table.loc[mask, column] = value
    save_frame(path, table, key='method', changed=set(display.values()), only_plinder=True)

tracking = PANEL / 'molecule-tracking'
labels = {mapping.get(m, m): label.replace('†', '') for m, label in json.loads((tracking / 'methods.json').read_text()).items()}
means = pd.read_csv(tracking / 'molecule_matching_values.csv').groupby('label', sort=False)[['mat_r', 'mat_p']].mean().sort_values('mat_r')
(tracking / 'methods.json').write_text(json.dumps({m: labels[m] for m in means.index}, indent=2) + '\n')
means.to_csv(tracking / 'matching_means.csv')

# Archive compact executed seed provenance, retaining hashes of complete records.
records = []
generation = benchmark / 'data/results/generation/plinder_independent_torsion'
for path in sorted((generation / 'provenance').glob('*.json')):
    info = json.loads(record(path).read_text())
    assert info['experimental_coordinates_used'] is False and info['status'] == 'ok'
    assert info['raw_conformer_count'] == info['minimized_conformer_count'] == 1000
    assert len(info['starts']) == 5 and all(s['retained'] == 200 for s in info['starts'])
    records.append({k: v for k, v in info.items() if k != 'mmff_statuses'})
assert len(records) == 22
(HERE / 'seed-provenance.json').write_text(json.dumps(records, indent=2) + '\n')
for relative in ['scripts/run_plinder_independent_torsion.py', 'src/casf_benchmark/generation/plinder_independent_torsion.py']:
    shutil.copy2(record(benchmark / relative), HERE / Path(relative).name)
audit = record(benchmark / 'data/mapping/casf_torsion_seed_manifest.csv')
audit_text = audit.read_text().replace('\\r\\n', '\n')
audit_rows = pd.read_csv(io.StringIO(audit_text))
assert audit_rows.groupby('ligand_set').size().to_dict() == {'core': 94, 'ref': 1236}
assert audit_rows.seed_source.eq('chembl3d_topology_sdf').all() and audit_rows.status.eq('unchanged').all()
(HERE / 'casf-seed-audit.csv').write_text(audit_text)
(HERE / 'provenance.json').write_text(json.dumps(dict(
    methods=display, replaced_methods=mapping, source_sha256=sources,
    scope='PLINDER torsion raw/minimized only; CASF and all other methods unchanged',
    protocol='5 ETKDGv3 starts from isomeric SMILES, 200 accepted torsion samples/start; no experimental coordinates',
    scoring='Published full-pool COV/MAT and clustering; PB separate; inclusive RMSD cutoffs',
    reference_counts=dict(molecules=22, conformations=1067), structural_recalculation=False), indent=2) + '\n')
print(summary[summary.label.isin(mapping.values())][['label','cov_r_0p75_mean','cov_p_0p75_mean','mat_r_mean','mat_p_mean']].to_string(index=False))
