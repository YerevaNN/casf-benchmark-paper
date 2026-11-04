"""Verify replacement measurements, unchanged comparators, and current table ranks."""
from pathlib import Path
import csv
import io
import json
import re
import subprocess

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
BASE = 'c0f04c4'
OLD = {'torsion_raw_druglike', 'torsion_minimized_22_druglike'}
NEW = {'torsion_independent_raw_druglike', 'torsion_independent_minimized_druglike'}

for relative, key in [
    ('plinder22-data/per-molecule.csv', 'label'),
    ('plinder22-data/summary.csv', 'label'),
    ('plinder22-data/threshold-curves.csv', 'method'),
    ('results-data/druglike_selected.csv', 'label'),
    ('plinder22-data/molecule-tracking/molecule_matching_values.csv', 'label'),
]:
    old = subprocess.check_output(['git', 'show', f'{BASE}:jcim-overleaf/{relative}'], cwd=PAPER, text=True)
    a = [r for r in csv.DictReader(io.StringIO(old)) if r[key] not in OLD]
    with (PAPER / relative).open() as stream:
        b = [r for r in csv.DictReader(stream) if r[key] not in NEW]
    assert a == b, relative
for number in [1, 2, 3]:
    relative = f'results-data/three-dataset-tables/table{number}.csv'
    old = subprocess.check_output(['git', 'show', f'{BASE}:jcim-overleaf/{relative}'], cwd=PAPER, text=True)
    a = [r for r in csv.DictReader(io.StringIO(old)) if r['cohort'] != 'plinder']
    with (PAPER / relative).open() as stream:
        b = [r for r in csv.DictReader(stream) if r['cohort'] != 'plinder']
    assert a == b, relative

current = pd.read_csv(PAPER / 'plinder22-data/per-molecule.csv')
for method in NEW:
    source = pd.read_csv(HERE / f'{method}-per_molecule.csv').set_index('smiles').sort_index()
    group = current[current.label.eq(method)].set_index('smiles').sort_index()
    pd.testing.assert_frame_equal(group[source.columns], source, check_dtype=False)
    assert len(group) == 22 and group.num_true_confs.sum() == 1067
    manifest = json.loads((HERE / f'{method}-generation.json').read_text())
    assert manifest['experimental_coordinates_used'] is False
    assert manifest['n_etkdg_starts'] == 5 and manifest['conformers_per_start'] == 200
    assert manifest['total_conformers'] == 22000
assert not current.label.isin(OLD).any()
assert all(r['experimental_coordinates_used'] is False for r in json.loads((HERE / 'seed-provenance.json').read_text()))

for number, filename, columns in [
    (1, 'diversity', ['clusters_0p5', 'clusters_1p0', 'largest_percent']),
    (2, 'count-recovery', ['recovery_0p5', 'recovery_0p75', 'distance']),
]:
    data = pd.read_csv(PAPER / f'results-data/three-dataset-tables/table{number}.csv')
    tex = (PAPER / f'tables/{filename}.tex').read_text()
    for ci, cohort in enumerate(['plinder', 'core', 'ref']):
        peers = data[data.cohort.eq(cohort) & data.method.ne('ChEMBL3D-PB')]
        for mi, metric in enumerate(columns):
            ranks = sorted(set(peers[metric].dropna()), reverse=metric not in ['distance', 'largest_percent'])[:2]
            for row in peers.to_dict('records'):
                line = next(l for l in tex.splitlines() if l.startswith(row['method'] + ' &'))
                cell = line.split(' & ')[1 + ci * 3 + mi]
                assert (r'\textbf{' in cell) == (row[metric] == ranks[0])
                assert (r'\underline{' in cell) == (len(ranks) > 1 and row[metric] == ranks[1])

means = pd.read_csv(PAPER / 'results-data/druglike_means.csv')
tex = (PAPER / 'tables/druglike.tex').read_text()
columns = ['cov_r_075', 'cov_p_075', 'mat_r', 'mat_p', 'num_gen_confs', 'pooled_pb_pass']
for index, metric in enumerate(columns):
    ranks = sorted(set(means[metric]), reverse=metric not in ['mat_r', 'mat_p'])[:2]
    for row in means.to_dict('records'):
        line = next(l for l in tex.splitlines() if l.startswith(row['display_label'] + ' &'))
        cell = line.split(' & ')[1 + index].replace(r' \\', '').strip()
        lower = metric in ['mat_r', 'mat_p']
        expected = f'{row[metric]:.3f}' if lower else f'{row[metric] * (1 if metric == "num_gen_confs" else 100):.1f}'
        assert re.sub(r'\\(?:textbf|underline)\{([^}]+)\}', r'\1', cell) == expected
        assert (r'\textbf{' in cell) == (metric != 'num_gen_confs' and np.isclose(row[metric], ranks[0], rtol=0, atol=1e-12))
        assert (r'\underline{' in cell) == (metric != 'num_gen_confs' and len(ranks) > 1 and np.isclose(row[metric], ranks[1], rtol=0, atol=1e-12))

for name in ['results.tex', 'methods.tex', 'discussion.tex', 'figures.tex',
             'tables/diversity.tex', 'tables/count-recovery.tex', 'tables/druglike.tex',
             'tables/supporting-druglike.tex', 'tables/generation-methods.tex']:
    text = (PAPER / name).read_text()
    assert 'dagger' not in text and '†' not in text and 'reference-seeded' not in text.lower(), name
print('Verified both independent torsion runs, unchanged other methods/CASF cells, all highlighted ranks, and removal of reference-seed markers.')
