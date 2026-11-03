"""Merge recalculated primary-Qwen energies and sampling into current paper inputs."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent;P=HERE.parents[1]
OLD='qwen_1p7b_fsq_bigdata_step47023';NEW='qwen_1p7b_v2_coordtoken_step41054'
def replace(path,frame,key='method'):
 with path.open() as f:r=csv.DictReader(f);fields=r.fieldnames;rows=list(r)
 kept=[r for r in rows if not r[key].startswith((OLD,NEW))]
 position=next(i for i,r in enumerate(rows) if r[key].startswith((OLD,NEW)))
 new=frame.copy()
 if 'label' in fields:new['label']='Qwen 1.7B FSQ'
 for col in fields:
  if col not in new:new[col]=''
 values=new[fields].fillna('').to_dict('records');kept[position:position]=values
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(kept)
for cohort in ['core','ref']:
 dest=P/'energy-data' if cohort=='core' else P/'results-data/core-ref-protocol/ref'
 e=pd.read_csv(HERE/cohort/'energy-per-molecule.csv')
 replace(dest/'per_molecule.csv',e)
 all_e=pd.read_csv(dest/'per_molecule.csv');b=all_e[all_e.method.eq('chembl3d_gt_pb')].set_index('mol_id')
 d=e[e.n_energy.gt(0)]
 summary=pd.DataFrame([dict(method=NEW+'_chembl_count',label='Qwen 1.7B FSQ',n_molecules=len(d),n_conformers=d.n_energy.sum(),n_failed=e.n_failed.sum(),mean_clusters_0p5=d.clusters_0p5.mean(),median_energy_std=d.std_energy.median(),mean_energy_std=d.std_energy.mean(),median_mean_energy=d.mean_energy.median(),median_mean_energy_delta_chembl=(d.set_index('mol_id').mean_energy-b.mean_energy).median(),median_archived_std=d.archived_std.median())])
 replace(dest/'summary.csv',summary)
 if cohort=='core':
  pc=pd.read_csv(dest/'per_conformer.csv.gz');new=pd.read_csv(HERE/cohort/'energy-per-conformer.csv.gz')
  pd.concat([pc[~pc.method.str.startswith((OLD,NEW))],new],ignore_index=True).to_csv(dest/'per_conformer.csv.gz',index=False)
 sampling=P/'sampling-data' if cohort=='core' else dest
 for name in ['recovery_curves.csv','recovery80_frontier.csv']:
  replace(sampling/name,pd.read_csv(HERE/cohort/name))
 if cohort=='core':
  replace(sampling/'endpoints.csv',pd.read_csv(HERE/cohort/'endpoints.csv'))
  for filename,arrays in [('candidate_rmsds.npz',['rmsd','n_candidates']),('entry_probabilities.npz',['probabilities'])]:
   z=np.load(sampling/filename);v={k:z[k].copy() for k in z.files};newz=np.load(HERE/cohort/filename)
   np.testing.assert_array_equal(v['mol_ids'],newz['mol_ids']);i=next(i for i,m in enumerate(v['methods']) if str(m).startswith((OLD,NEW)))
   # Widen string dtype before substituting a longer checkpoint identity.
   methods=[str(m) for m in v['methods']];methods[i]=NEW+'_fixed';v['methods']=np.array(methods)
   for a in arrays:v[a][i]=newz[a][0]
   np.savez_compressed(sampling/filename,**v)
 else:
  c=pd.read_csv(HERE/cohort/'recovery_curves.csv');f=pd.read_csv(HERE/cohort/'recovery80_frontier.csv')
  end=pd.DataFrame([dict(method=NEW+'_fixed',recovery=100*c[c.budget.eq(1000)&c.threshold.eq(.75)].recovery.iloc[0],rmsd_threshold=f[f.budget.eq(1000)].rmsd_threshold.iloc[0])])
  replace(dest/'full-pool-endpoints.csv',end)
 print(cohort,summary.to_string(index=False));print('energy totals',int(all_e.n_energy.sum()),int(all_e.n_retained.sum()))
# Update ref protocol counts and provenance; retain original calculation lineage.
p=P/'results-data/ref-panel/counts.json';counts=json.loads(p.read_text());e=pd.read_csv(P/'results-data/core-ref-protocol/ref/per_molecule.csv')
counts.update(retained_energy_conformers=int(e.n_retained.sum()),finite_energy_conformers=int(e.n_energy.sum()),undefined_energies=int(e.n_failed.sum()),empty_energy_ensembles=int(e.n_energy.eq(0).sum()))
p.write_text(json.dumps(counts,indent=2)+'\n')
for rel in ['results-data/ref-panel/provenance.json','results-data/core-ref-protocol/ref/provenance.json','energy-data/provenance.json','sampling-data/extraction-provenance.json','plinder22-data/provenance.json','plinder22-data/molecule-tracking/source-provenance.json']:
 p=P/rel;v=json.loads(p.read_text());v['primary_qwen_update']={'checkpoint':NEW,'provenance':'results-data/qwen-v2/provenance.json','calculation_provenance':'results-data/qwen-v2/calculation-provenance.json'}
 if rel=='results-data/ref-panel/provenance.json':v['counts']=counts;v['methods']['Qwen 1.7B FSQ']=NEW;v['core_and_plinder_preserved']='At original ref update; primary Qwen subsequently replaced in all three panels.'
 if rel=='plinder22-data/provenance.json':v['methods']=json.loads((P/'plinder22-data/methods.json').read_text());v['qwen']=NEW
 p.write_text(json.dumps(v,indent=2)+'\n')
