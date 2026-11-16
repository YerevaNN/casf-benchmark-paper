"""Check the archived Figure 4 curves, full-pool endpoints and rendered artifacts."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
provenance = json.loads((HERE/'provenance.json').read_text())
for name, digest in provenance['archive_sha256'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest, name
rendered = json.loads((HERE/'rendered-files.json').read_text())
for section in ['inputs', 'outputs']:
    for name, digest in rendered[section].items():
        assert hashlib.sha256((PAPER/name).read_bytes()).hexdigest() == digest, name

per = pd.read_csv(HERE/'plinder-per-molecule-curves.csv')
summary = pd.read_csv(HERE/'plinder-curves.csv')
hits = pd.read_csv(HERE/'plinder-reference-hit-counts.csv')
published = pd.read_csv(PAPER/'plinder22-data/per-molecule.csv')
casf = pd.read_csv(PAPER/'sampling-data/recovery_curves.csv')
budgets = sorted(casf.budget.unique())
assert len(budgets) == 63 and sorted(per.budget.unique()) == budgets
assert len(per) == 4*22*63 and len(hits) == 4*1067
assert not per.duplicated(['method','smiles','budget']).any()
assert not hits.duplicated(['method','smiles','reference_index']).any()
assert per.method.nunique() == 4
assert 'qwen_1p7b_v2_coordtoken_step41054' in set(per.method)
np.testing.assert_array_equal(per.n_sampled, np.minimum(per.budget, per.n_generated))
for (method, smiles), rows in per.groupby(['method','smiles']):
    rows = rows.sort_values('budget')
    assert (np.diff(rows.coverage) >= -1e-12).all()
    assert (np.diff(rows.mat_r) <= 1e-12).all()
    target = published[published.label.eq(method) & published.smiles.eq(smiles)].iloc[0]
    assert rows.n_generated.nunique() == 1
    assert rows.n_generated.iloc[0] == target.n_generated
    np.testing.assert_allclose(rows.iloc[-1][['coverage','mat_r']].to_numpy(dtype=float),
                               [target.cov_r_0p75,target.mat_r],atol=1e-12,rtol=0)
    h = hits[hits.method.eq(method) & hits.smiles.eq(smiles)]
    assert len(h) == target.num_true_confs
    assert (h.n_generated == target.n_generated).all()
    expected = [hypergeom.sf(0,int(target.n_generated),h.n_hits.to_numpy(),int(k)).mean()
                for k in rows.n_sampled]
    np.testing.assert_allclose(rows.coverage,expected,atol=1e-12,rtol=0)
means = per.groupby(['method','budget'])[['coverage','mat_r']].mean().sort_index()
np.testing.assert_allclose(means,summary.set_index(['method','budget']).sort_index()[means.columns],atol=1e-12,rtol=0)
assert (summary.coverage.between(0,1)).all() and (summary.mat_r>=0).all()
print('Verified 88 full-pool COV-R/MAT-R pairs, 4,268 reference hit counts, all 63 budgets, monotonicity, weighting, and figure hashes.')
