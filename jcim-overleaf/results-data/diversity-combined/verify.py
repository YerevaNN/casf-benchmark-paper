"""Verify Figure 1/S8 values, artifact hashes and current manuscript inclusion."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
PAPER=HERE.parents[1]
provenance=json.loads((HERE/'provenance.json').read_text())
for section in ['inputs','outputs']:
 for name,digest in provenance[section].items():
  assert hashlib.sha256((PAPER/name).read_bytes()).hexdigest()==digest,name
pl=pd.read_csv(HERE/'plinder-points.csv').set_index('method')
original_pl=pd.read_csv(PAPER/'plinder22-data/scatter-points.csv').set_index('method')
assert len(pl)==10 and set(pl.index)==set(original_pl.index)
np.testing.assert_allclose(pl[['clusters_0p5','mat_r']],original_pl.loc[pl.index,['clusters_0p5','mat_r']],atol=1e-12,rtol=0)
casf=pd.read_csv(HERE/'casf-points.csv').set_index(['cohort','method'])
original_casf=pd.read_csv(PAPER/'results-data/three-dataset-tables/casf-scatter-points.csv').set_index(['cohort','method'])
assert len(casf)==22 and set(casf.index)==set(original_casf.index)
np.testing.assert_allclose(casf[['clusters_0p5','distance']],original_casf.loc[casf.index,['clusters_0p5','distance']],atol=1e-12,rtol=0)
curves=pd.read_csv(PAPER/'plinder22-data/threshold-curves.csv')
summary=pd.read_csv(PAPER/'plinder22-data/summary.csv').set_index('label')
for method in ['rdkit_random_raw_druglike','mcf_drugs_l_druglike','flowr_druglike','torsional_diffusion_druglike','nextmol_dmt_l_druglike','loqi_druglike','qwen_1p7b_v2_coordtoken_step41054']:
 row=curves[curves.method.eq(method)&curves.threshold.eq(.75)].iloc[0]
 # The threshold arrays must remain unmodified; the renderer reads them directly.
 assert np.isfinite([row.recall,row.precision]).all()
 assert (curves[curves.method.eq(method)][['recall','precision']].diff().dropna()>=-1e-12).all().all()
results=(PAPER/'results.tex').read_text()
assert r'\drugfigure' not in results and 'fig:druglike' not in results
assert results.count(r'\diversityfigure')==1
assert 'fig:plinder-thresholds' in results
assert 'plinder-thresholds-appendix' in (PAPER/'supporting-content.tex').read_text()
print('Verified all 32 scatter coordinates unchanged, seven threshold curves, figure hashes and manuscript inclusion.')
