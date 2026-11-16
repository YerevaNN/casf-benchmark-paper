"""Validate the archived step-64000 measurements, unused in current comparisons."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
PAPER=HERE.parents[1]
NEW='qwen_0p6b_casf_dedup_binned_pm12_step64000'
d=pd.read_csv(HERE/'per-molecule.csv').set_index('smiles')
assert len(d)==22 and d.index.is_unique and d.num_true_confs.sum()==1067
with np.load(HERE/'nearest-distances.npz',allow_pickle=False) as archive:
 for entry in json.loads((HERE/'entries.json').read_text()):
  row=d.loc[entry['smiles']];key=entry['key']
  r,g=archive[key+'_reference'],archive[key+'_generated']
  assert len(r)==row.num_true_confs and len(g)==row.num_gen_confs
  assert np.isfinite(r).all() and np.isfinite(g).all()
  np.testing.assert_allclose([r.mean(),g.mean()],row[['mat_r','mat_p']].astype(float),atol=2e-6,rtol=0)
  np.testing.assert_allclose([(r<=.75).mean(),(g<=.75).mean()],row[['cov_r_0p75','cov_p_0p75']].astype(float),atol=1e-12,rtol=0)
print('Verified the archived step-64000 scores; this run is excluded from current manuscript comparisons.')
