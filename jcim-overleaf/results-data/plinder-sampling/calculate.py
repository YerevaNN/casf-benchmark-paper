"""Recalculate PLINDER curves using the original RMSD caches listed in provenance.json."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from scipy.stats import hypergeom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from casf_benchmark.analysis.druglike_covmat import load_rmsd_matrices
OUT=Path(__file__).resolve().parent
P=OUT.parents[1]
STYLE=json.loads((P/'sampling-data/palette.json').read_text())
methods=json.loads((P/'sampling-data/display.json').read_text())['methods']
pl_methods=['qwen_1p7b_v2_coordtoken_step41054','loqi_druglike','flowr_druglike','nextmol_dmt_l_druglike']
names=['Qwen','LoQI','FlowR','NExT-Mol']
colors=[STYLE[m.removesuffix('_fixed')][1] for m in methods]
core=pd.read_csv(P/'sampling-data/recovery_curves.csv')
ref=pd.read_csv(P/'results-data/core-ref-protocol/ref/recovery_curves.csv')
budgets=sorted(core.budget.unique()); assert budgets[0]==10 and budgets[-1]==1000
per=pd.read_csv(P/'plinder22-data/per-molecule.csv')
provenance=json.loads((OUT/'provenance.json').read_text())['sources']
sources={}; rows=[]; hit_counts=[]

def expected_minimum(sorted_distances, k):
 # Probability that rank j (zero-based) is the smallest of k sampled columns.
 n=sorted_distances.shape[1]; k=min(k,n)
 weights=np.zeros(n); weights[0]=k/n
 for j in range(1,n-k+1):
  weights[j]=weights[j-1]*(n-j-k+1)/(n-j)
 assert np.isclose(weights.sum(),1,atol=1e-12)
 support=n-k+1
 assert np.isfinite(sorted_distances[:,:support]).all()
 return float((sorted_distances[:,:support] @ weights[:support]).mean())

for method in pl_methods:
 path=Path(next(path for path in provenance if path.endswith('/'+method+'.npz')))
 digest=hashlib.sha256(path.read_bytes()).hexdigest();assert digest==provenance[str(path)]
 sources[str(path)]=digest
 matrices=load_rmsd_matrices(path); current=per[per.label.eq(method)].set_index('smiles')
 assert set(current.index)==set(matrices) and len(matrices)==22
 for smiles,matrix in matrices.items():
  n=matrix.shape[1]; assert matrix.shape==(current.loc[smiles,'num_true_confs'],current.loc[smiles,'n_generated'])
  assert (~np.isfinite(matrix)).sum(axis=1).max()<min(budgets)
  ordered=np.sort(np.where(np.isfinite(matrix),matrix,np.inf),axis=1)
  if np.isfinite(matrix).all():
   assert np.isclose(expected_minimum(ordered,1),matrix.mean())
  assert np.isclose(expected_minimum(ordered,n),current.loc[smiles,'mat_r'])
  hits=(matrix<=.75).sum(axis=1)
  full=(hits>0).mean();assert abs(full-current.loc[smiles,'cov_r_0p75'])<1e-12
  for ri,h in enumerate(hits):hit_counts.append(dict(method=method,smiles=smiles,reference_index=ri,n_generated=n,n_hits=int(h)))
  for k in budgets:
   prob=hypergeom.sf(0,n,hits,min(k,n))
   rows.append(dict(method=method,smiles=smiles,name=current.loc[smiles,'name'],budget=k,n_generated=n,n_sampled=min(k,n),coverage=float(prob.mean()),mat_r=expected_minimum(ordered,k)))
d=pd.DataFrame(rows);d.to_csv(OUT/'plinder-per-molecule-curves.csv',index=False)
pd.DataFrame(hit_counts).to_csv(OUT/'plinder-reference-hit-counts.csv',index=False)
summary=d.groupby(['method','budget'],as_index=False)[['coverage','mat_r']].mean();summary.to_csv(OUT/'plinder-curves.csv',index=False)
for method,g in d.groupby('method'):
 table=g.pivot(index='smiles',columns='budget',values='coverage')
 assert (np.diff(table.values,axis=1)>=-1e-12).all()
 assert np.allclose(table[1000].values,per[per.label.eq(method)].set_index('smiles').loc[table.index,'cov_r_0p75'])
for method,g in d.groupby('method'):
 table=g.pivot(index='smiles',columns='budget',values='mat_r')
 assert (np.diff(table.values,axis=1)<=1e-12).all()
 assert np.allclose(table[1000].values,per[per.label.eq(method)].set_index('smiles').loc[table.index,'mat_r'])
print("Verified and recalculated PLINDER coverage and MAT-R curves from archived RMSD matrices.")
