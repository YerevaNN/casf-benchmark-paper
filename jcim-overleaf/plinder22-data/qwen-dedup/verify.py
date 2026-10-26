"""Validate the additional Qwen row and the approved molecule-tracking figure."""
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
means=pd.read_csv(PAPER/'results-data/druglike_means.csv').set_index('label')
assert len(means)==11
for dest,src in [('mat_r','mat_r'),('mat_p','mat_p'),('cov_r_075','cov_r_0p75'),('cov_p_075','cov_p_0p75'),('num_gen_confs','num_gen_confs')]:
 np.testing.assert_allclose(means.loc[NEW,dest],d[src].mean(),atol=1e-12,rtol=0)
np.testing.assert_allclose(means.loc[NEW,'pooled_pb_pass'],d.n_pb_pass.sum()/d.n_generated.sum(),atol=1e-12,rtol=0)
tex=(PAPER/'tables/druglike.tex').read_text()
line=next(l for l in tex.splitlines() if l.startswith('Qwen 0.6B CASF-dedup &'))
assert line==r'Qwen 0.6B CASF-dedup & 89.3 & 31.8 & 0.441 & 1.113 & 997.8 & 60.5 \\'
F=PAPER/'plinder22-data/molecule-tracking'
plot=pd.read_csv(F/'molecule_matching_values.csv')
base=pd.read_csv(PAPER/'plinder22-data/per-molecule.csv')
allrows=pd.concat([base,d.reset_index()],ignore_index=True)
matched=plot.merge(allrows,on=['label','smiles'],validate='one_to_one',suffixes=('_plot','_source'))
assert len(matched)==198 and plot.smiles.nunique()==22
for col in ['mat_r','mat_p','num_true_confs','num_gen_confs']:
 np.testing.assert_allclose(matched[col+'_plot'],matched[col+'_source'],atol=1e-12,rtol=0)
for name,sha in json.loads((F/'paper-artifacts.json').read_text()).items():
 assert hashlib.sha256((PAPER/name).read_bytes()).hexdigest()==sha,name
print('Verified the new Qwen scores, pooled validity, Table 4 row, 198 plotted records, and figure hashes.')
