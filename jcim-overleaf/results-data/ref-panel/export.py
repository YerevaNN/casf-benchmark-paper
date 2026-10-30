"""Refresh manuscript ref inputs from the published 1,044-crystal cohort.

No geometry, energy, generation, or sampling calculations are rerun.
"""
from pathlib import Path
import argparse,csv,hashlib,json,shutil,sqlite3
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
PAPER=HERE.parents[1]
BENCH=Path('/mnt/weka/mbedrosian/code/casf-benchmark')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--database',type=Path,default=BENCH/'data/results/casf_analysis_dashboard.sqlite')
parser.add_argument('--manifest',type=Path,default=BENCH/'data/mapping/casf16_ref_dedup_selection_v1.csv')
args=parser.parse_args()
manifest=pd.read_csv(args.manifest)
selected=manifest.loc[manifest.selected.astype(str).str.lower().eq('true')]
ids=set(selected.ligand_id)
assert len(ids)==len(selected)==selected.casf_heavy_isomeric_smiles.nunique()==1044
styles=json.loads((PAPER/'sampling-data/palette.json').read_text())
name_to_key={v[0]:k for k,v in styles.items()}
T=PAPER/'results-data/three-dataset-tables'
names=list(pd.read_csv(T/'table1.csv').method.unique())
keys=[name_to_key[n] for n in names]
methods=[k+'_'+b for k in keys if k!='chembl3d_gt_pb' for b in ['chembl_count','fixed']]+['chembl3d_gt_pb']
columns=['method','mol_id','heavy_atoms','rotatable_bonds','casf_best_rmsd','casf_hit_0p5','casf_hit_0p75','greedy_clusters_0p5','greedy_clusters_1p0','largest_cluster_fraction_1p0','post_pb_confs']
with sqlite3.connect(f'file:{args.database}?mode=ro',uri=True) as c:
 per=pd.read_sql_query('SELECT '+','.join(columns)+" FROM per_ligand_long WHERE ligand_set='ref_dedup'",c)
 per=per[per.method.isin(methods)].sort_values(['method','mol_id'])
 assert len(per)==1044*21 and not per.duplicated(['method','mol_id']).any()
 for method,d in per.groupby('method'):
  assert set(d.mol_id)==ids
  np.testing.assert_allclose(d.casf_best_rmsd.le(.5).astype(float),d.casf_hit_0p5.fillna(0),atol=0)
  np.testing.assert_allclose(d.casf_best_rmsd.le(.75).astype(float),d.casf_hit_0p75.fillna(0),atol=0)
 def read(table):
  return pd.read_sql_query(f"SELECT * FROM {table} WHERE ligand_set='ref_dedup'",c).drop(columns='ligand_set')
 energy=read('paper_energy_summary');energy_per=read('paper_energy_per_ligand')
 curves=read('paper_sampling_curves');frontier=read('paper_sampling_frontier');reference=read('paper_sampling_reference')
 protocol=read('paper_protocol_runs')
per.to_csv(HERE/'per-entry.csv',index=False)
shutil.copy2(args.manifest,HERE/'selection.csv')
selection_info=json.loads(args.manifest.with_suffix('.json').read_text())
(HERE/'selection.json').write_text(json.dumps(selection_info,indent=2)+'\n')
coremap=pd.read_csv(BENCH/'data/mapping/casf16_core_chembl3d_exact_intersection.csv')
# Reuse the stored numerical strings for core and PLINDER, replacing ref rows only.
changes=[]
for number in [1,2,3]:
 path=T/f'table{number}.csv'
 with path.open() as f: reader=csv.DictReader(f);fields=reader.fieldnames;rows=list(reader)
 for row in rows:
  if row['cohort']!='ref':continue
  key=name_to_key[row['method']]
  method=key if key=='chembl3d_gt_pb' else key+'_'+(row['budget'] if number==3 else 'chembl_count')
  d=per[per.method.eq(method)];assert len(d)==1044
  old=dict(row)
  row['n_total']='1044'
  if number==1:
   row.update(n_defined=str(d.greedy_clusters_0p5.count()),clusters_0p5=str(d.greedy_clusters_0p5.mean()),clusters_1p0=str(d.greedy_clusters_1p0.mean()),largest_percent=str(100*d.largest_cluster_fraction_1p0.mean()),status='published_ref_selection')
  elif number==2:
   row.update(n_distance=str(d.casf_best_rmsd.count()),recovery_0p5=str(100*d.casf_best_rmsd.le(.5).mean()),recovery_0p75=str(100*d.casf_best_rmsd.le(.75).mean()),distance=str(d.casf_best_rmsd.mean()))
  else:
   row.update(n_defined=str(d.casf_best_rmsd.count()),distance=str(d.casf_best_rmsd.mean()),recovery=str(100*d.casf_best_rmsd.le(.75).mean()),clusters=str(d.greedy_clusters_0p5.mean()),mean_retained=str(d.post_pb_confs.fillna(0).mean()),status='published_ref_selection')
  for metric in row:
   if metric in ['n_total','n_defined','n_distance','clusters_0p5','clusters_1p0','largest_percent','recovery_0p5','recovery_0p75','distance','recovery','clusters','mean_retained']:
    changes.append(dict(table=number,method=row['method'],budget=row.get('budget','chembl_count'),metric=metric,old=float(old[metric]),new=float(row[metric]),delta=float(row[metric])-float(old[metric])))
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
if not (HERE/'changes-from-1236.csv').exists():
 pd.DataFrame(changes).to_csv(HERE/'changes-from-1236.csv',index=False)
# Full-protocol energy and sampling tables are already recomputed upstream on this selection.
E=PAPER/'results-data/core-ref-protocol/ref'
energy.to_csv(E/'summary.csv',index=False);energy_per.to_csv(E/'per_molecule.csv',index=False)
curves.to_csv(E/'recovery_curves.csv',index=False);frontier.to_csv(E/'recovery80_frontier.csv',index=False)
refdict={'hit_rates':{str(r.threshold):r.recovery for r in reference.itertuples()},'target_recovery':float(reference.target_recovery.iloc[0]),'rmsd_threshold':float(reference.rmsd_threshold.iloc[0]),'n_entries':1044}
(E/'reference.json').write_text(json.dumps(refdict,indent=2)+'\n')
endmethods=['qwen_1p7b_fsq_bigdata_step47023_fixed','loqi_raw_fixed','flowr_raw_fixed','nextmol_dmt_l_raw_fixed']
ends=[]
for m in endmethods:
 d=per[per.method.eq(m)]
 end=curves[curves.method.eq(m)&curves.budget.eq(1000)&curves.threshold.eq(.75)].iloc[0]
 assert abs(end.recovery-d.casf_best_rmsd.le(.75).mean())<1e-12
 threshold=frontier[frontier.method.eq(m)&frontier.budget.eq(1000)].rmsd_threshold.iloc[0]
 ends.append(dict(method=m,recovery=100*end.recovery,rmsd_threshold=threshold))
ends.append(dict(method='chembl3d_gt_pb',recovery=100*refdict['hit_rates']['0.75'],rmsd_threshold=refdict['rmsd_threshold']))
pd.DataFrame(ends).to_csv(E/'full-pool-endpoints.csv',index=False)
base=per[per.method.eq('chembl3d_gt_pb')].set_index('mol_id')
S=PAPER/'results-data/size-flexibility'
for filename,descriptor,labels,bins in [('casf-size-source.csv','heavy_atoms',['<20','20-29','30-39','40+'],[0,20,30,40,np.inf]),('casf-flexibility-source.csv','rotors',['0-3','4-6','7-8','9+'],[-1,3,6,8,np.inf])]:
 path=S/filename
 with path.open() as f:reader=csv.DictReader(f);fields=reader.fieldnames;rows=[r for r in reader if r['set']!='ref']
 grouping=pd.cut(base.heavy_atoms if descriptor=='heavy_atoms' else base.rotatable_bonds,bins,labels=labels,right=descriptor!='heavy_atoms')
 assert grouping.notna().all()
 for key in keys:
  m=key if key=='chembl3d_gt_pb' else key+'_fixed';d=per[per.method.eq(m)].set_index('mol_id').reindex(base.index)
  for label in labels:
   g=d[grouping.eq(label)];n=len(g);hits=int(g.casf_best_rmsd.le(.75).sum())
   row=dict(set='ref',method=m,n=n,hit_075=hits/n if n else '',**{descriptor:label})
   if 'hits' in fields:row['hits']=hits
   if 'rmsd' in fields:row['rmsd']=g.casf_best_rmsd.mean()
   rows.append(row)
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
summary={'entries':1044,'unique_stereochemical_smiles':1044,'parent_compounds':int(selected.chembl3d_mol_id.nunique()),'shared_parent_compounds_core':len(set(selected.chembl3d_mol_id)&set(coremap.chembl3d_mol_id)),'measurable_count_min':int(per[per.method.str.endswith('_chembl_count')].groupby('method').casf_best_rmsd.count().min()),'measurable_count_max':1044,'retained_energy_conformers':int(energy_per.n_retained.sum()),'finite_energy_conformers':int(energy_per.n_energy.sum()),'undefined_energies':int(energy_per.n_failed.sum()),'empty_energy_ensembles':int(energy_per.n_energy.eq(0).sum())}
(HERE/'counts.json').write_text(json.dumps(summary,indent=2)+'\n')
prov={'database':str(args.database.resolve()),'database_sha256':hashlib.sha256(args.database.read_bytes()).hexdigest(),'database_cohort':'ref_dedup','manuscript_cohort':'ref','selection':selection_info,'counts':summary,'protocol':protocol.to_dict('records'),'methods':dict(zip(names,keys)),'geometric_recalculation':False,'core_and_plinder_preserved':True}
(HERE/'provenance.json').write_text(json.dumps(prov,indent=2)+'\n')
oldprov=json.loads((E/'provenance.json').read_text())
oldprov['ref_selection']={'source':'results-data/ref-panel/provenance.json','n_entries':1044,'source_cohort':'ref_dedup','note':'Current ref is the frozen one-crystal-per-stereochemical-identity selection; original calculation lineage remains below.'}
(E/'provenance.json').write_text(json.dumps(oldprov,indent=2)+'\n')
print(json.dumps(summary,indent=2));print(pd.DataFrame(ends).to_string(index=False))
