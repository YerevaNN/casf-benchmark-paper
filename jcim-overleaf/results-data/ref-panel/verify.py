"""Verify current ref aggregates, source selection, and preserved other panels."""
from pathlib import Path
import csv,hashlib,io,json,re,subprocess,sys
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
PAPER=HERE.parents[1];REPO=PAPER.parent
BASELINE='5123902'
manifest=pd.read_csv(HERE/'selection.csv')
selected=manifest[manifest.selected.astype(str).str.lower().eq('true')]
assert len(selected)==selected.casf_heavy_isomeric_smiles.nunique()==1044
for row in manifest.itertuples():
 digest=hashlib.sha256(f'{row.selection_seed}\t{row.casf_heavy_isomeric_smiles}\t{row.ligand_id}'.encode()).hexdigest()
 assert digest==row.selection_digest
first=manifest.sort_values(['selection_digest','ligand_id']).groupby('casf_heavy_isomeric_smiles').first()
assert set(first.ligand_id)==set(selected.ligand_id)
per=pd.read_csv(HERE/'per-entry.csv')
assert len(per)==21924 and not per.duplicated(['method','mol_id']).any()
for _,d in per.groupby('method'):assert set(d.mol_id)==set(selected.ligand_id)
keys=json.loads((HERE/'provenance.json').read_text())['methods']
T=PAPER/'results-data/three-dataset-tables'
for number in [1,2,3]:
 path=T/f'table{number}.csv';data=pd.read_csv(path)
 old=list(csv.DictReader(io.StringIO(subprocess.check_output(['git','show',f'{BASELINE}:{path.relative_to(REPO)}'],cwd=REPO,text=True))))
 new=list(csv.DictReader(path.open()))
 assert [r for r in old if r['cohort']!='ref']==[r for r in new if r['cohort']!='ref']
 for row in data[data.cohort.eq('ref')].itertuples():
  key=keys[row.method];m=key if key=='chembl3d_gt_pb' else key+'_'+(row.budget if number==3 else 'chembl_count')
  d=per[per.method.eq(m)];assert row.n_total==1044
  if number==1:
   actual=[row.n_defined,row.clusters_0p5,row.clusters_1p0,row.largest_percent]
   expected=[d.greedy_clusters_0p5.count(),d.greedy_clusters_0p5.mean(),d.greedy_clusters_1p0.mean(),100*d.largest_cluster_fraction_1p0.mean()]
  elif number==2:
   actual=[row.n_distance,row.recovery_0p5,row.recovery_0p75,row.distance]
   expected=[d.casf_best_rmsd.count(),100*d.casf_best_rmsd.le(.5).mean(),100*d.casf_best_rmsd.le(.75).mean(),d.casf_best_rmsd.mean()]
  else:
   actual=[row.n_defined,row.distance,row.recovery,row.clusters,row.mean_retained]
   expected=[d.casf_best_rmsd.count(),d.casf_best_rmsd.mean(),100*d.casf_best_rmsd.le(.75).mean(),d.greedy_clusters_0p5.mean(),d.post_pb_confs.fillna(0).mean()]
  np.testing.assert_allclose(actual,expected,rtol=0,atol=1e-10)
E=PAPER/'results-data/core-ref-protocol/ref'
e=pd.read_csv(E/'per_molecule.csv');s=pd.read_csv(E/'summary.csv').set_index('method')
assert e.groupby('method').size().eq(1044).all()
base=e[e.method.eq('chembl3d_gt_pb')].set_index('mol_id').mean_energy
for m,d in e.groupby('method'):
 valid=d[d.n_energy.gt(0)]
 assert set(d.mol_id)==set(selected.ligand_id)
 row=s.loc[m]
 np.testing.assert_allclose([row.n_molecules,row.n_conformers,row.n_failed,row.median_energy_std,row.mean_energy_std,row.median_mean_energy,row.median_mean_energy_delta_chembl,row.mean_clusters_0p5],
 [len(valid),valid.n_energy.sum(),d.n_failed.sum(),valid.std_energy.median(),valid.std_energy.mean(),valid.mean_energy.median(),(valid.set_index('mol_id').mean_energy-base).median(),valid.clusters_0p5.mean()],rtol=1e-12,atol=1e-9)
assert e.n_energy.sum()==734206 and e.n_retained.sum()==734207 and e.n_failed.sum()==1
curves=pd.read_csv(E/'recovery_curves.csv');front=pd.read_csv(E/'recovery80_frontier.csv')
assert curves.recovery.between(0,1).all()
for (m,t),d in curves.groupby(['method','threshold']):
 assert (np.diff(d.sort_values('budget').recovery)>=-1e-12).all()
 expected=per[per.method.eq(m)].casf_best_rmsd.le(t).mean()
 np.testing.assert_allclose(d[d.budget.eq(1000)].recovery,expected,rtol=0,atol=1e-12)
for m,d in front.groupby('method'):
 finite=d.sort_values('budget').rmsd_threshold.dropna().to_numpy()
 assert (np.diff(finite)<=1e-12).all()
reference=json.loads((E/'reference.json').read_text());assert reference['n_entries']==1044
b=per[per.method.eq('chembl3d_gt_pb')].set_index('mol_id')
for t,val in reference['hit_rates'].items():np.testing.assert_allclose(val,b.casf_best_rmsd.le(float(t)).mean(),rtol=0,atol=1e-12)
S=PAPER/'results-data/size-flexibility'
for file,col,descriptor,bins,labels,right in [('casf-size-source.csv','heavy_atoms','heavy_atoms',[0,20,30,40,np.inf],['<20','20-29','30-39','40+'],False),('casf-flexibility-source.csv','rotors','rotatable_bonds',[-1,3,6,8,np.inf],['0-3','4-6','7-8','9+'],True)]:
 data=pd.read_csv(S/file);groups=pd.cut(b[descriptor],bins,labels=labels,right=right)
 for m,d in data[data['set'].eq('ref')].groupby('method'):
  assert d.n.sum()==1044
  scores=per[per.method.eq(m)].set_index('mol_id').reindex(b.index)
  for row in d.to_dict('records'):
   sub=scores[groups.eq(row[col])];assert row['n']==len(sub)
   np.testing.assert_allclose(row['hit_075'],sub.casf_best_rmsd.le(.75).mean(),rtol=0,atol=1e-12)
 # Core input strings remain exactly as before.
 old=list(csv.DictReader(io.StringIO(subprocess.check_output(['git','show',f'{BASELINE}:{(S/file).relative_to(REPO)}'],cwd=REPO,text=True))))
 new=list(csv.DictReader((S/file).open()))
 assert [r for r in old if r['set']=='core']==[r for r in new if r['set']=='core']
# No stale denominators or dashboard-only name in active manuscript source.
for path in PAPER.rglob('*.tex'):
 if 'preview' in path.parts:continue
 text=re.sub(r'(?<!\\)%[^\n]*','',path.read_text())
 assert not re.search(r'1,236|1236|1,220|835,97[23]|ref[ _-]dedup',text),(path,'stale ref description')
# Table 4 and PLINDER artwork remain byte-for-byte unchanged.
for name in ['tables/druglike.tex','tables/supporting-druglike.tex','results-data/druglike_means.csv','plinder22-data/per-molecule.csv','figures/figure-plinder-combined.pdf','figures/figure-plinder-molecule-tracking.pdf']:
 old=subprocess.check_output(['git','show',f'{BASELINE}:jcim-overleaf/{name}'],cwd=REPO)
 assert old==(PAPER/name).read_bytes(),name
# Table regeneration is deterministic and does not modify the reviewed source.
for script,names in [('results-data/three-dataset-tables/build_tables.py',['diversity','count-recovery','recovery']),('results-data/ref-panel/build_supporting.py',['supporting-recovery','size']),('results-data/core-ref-protocol/build_energy_table.py',['energy-summary'])]:
 before={n:(PAPER/f'tables/{n}.tex').read_bytes() for n in names}
 subprocess.run([sys.executable,str(PAPER/script)],check=True,capture_output=True)
 assert all((PAPER/f'tables/{n}.tex').read_bytes()==v for n,v in before.items()),script
print('Verified 1,044 unique selected identities, 21 ref pools, all ref table metrics, energy summaries, sampling endpoints, and size/flexibility strata. Core and PLINDER preserved; no stale ref counts in TeX.')
