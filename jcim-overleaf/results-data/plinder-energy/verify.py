"""Check Figure 2 against archived molecule summaries and rendering inputs."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
rows = pd.read_csv(HERE / 'per-molecule.csv')
source = pd.read_csv(PAPER / 'plinder22-data/per-molecule.csv')
keys = ['label', 'smiles']
pd.testing.assert_frame_equal(
    rows.set_index(keys).sort_index(),
    source[rows.columns].set_index(keys).sort_index(),
    check_dtype=False,
)
assert len(rows) == 220 and rows.label.nunique() == 10
np.testing.assert_array_equal(rows.energy_count, rows.n_pb_pass)
assert rows.energy_count.sum() == 192356
empty = rows[rows.energy_count.eq(0)]
assert empty[['label', 'name']].values.tolist() == [['flowr_druglike', '4-Hydroxytamoxifen']]
assert rows.energy_std.isna().equals(rows.energy_count.eq(0))
summary = pd.read_csv(HERE / 'summary.csv')
for r in summary.itertuples():
    group = rows[rows.label.eq(r.method)]
    measured = group[group.energy_count.gt(0)]
    assert len(group) == r.n_total == 22
    assert len(measured) == r.n_molecules
    np.testing.assert_allclose(
        [measured.energy_std.median(), measured.energy_std.mean(), measured.greedy_clusters_0p5.mean()],
        [r.median_energy_std, r.mean_energy_std, r.mean_clusters_0p5],
        rtol=1e-12,
    )
points = pd.read_csv(HERE / 'figure-points.csv')
assert len(points) == 64
for cohort, data in [
    ('plinder', summary),
    ('core', pd.read_csv(PAPER / 'energy-data/summary.csv')),
    ('ref', pd.read_csv(PAPER / 'results-data/core-ref-protocol/ref/summary.csv')),
]:
    for r in points[points.cohort.eq(cohort)].itertuples():
        original = data[data.label.eq(r.method)].iloc[0]
        np.testing.assert_allclose(r.value, original[r.metric], rtol=1e-12)
        assert r.n_molecules == original.n_molecules
manifest = json.loads((HERE / 'rendered-files.json').read_text())
for section in ['inputs', 'outputs']:
    for filename, digest in manifest[section].items():
        assert hashlib.sha256((PAPER / filename).read_bytes()).hexdigest() == digest, filename
audit = json.loads((HERE / 'protocol-audit.json').read_text())
assert hashlib.sha256((HERE / 'source-plinder-panel.py').read_bytes()).hexdigest() == audit['original_evaluator_sha256']
assert audit['protocol']['coordinate_optimization'] is False
for function in ['common_hydrogen_mmff_energy', 'energy_summary']:
    assert all(audit['function_comparisons'][function].values())
print('Verified 220 molecule records, 192,356 energies, 64 plotted values and rendering hashes.')
