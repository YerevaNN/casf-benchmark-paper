"""Replace the primary FSQ checkpoint in active paper inputs; retain other methods."""
from pathlib import Path
import csv,hashlib,json,shutil,sqlite3
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent;P=HERE.parents[1]
B=Path('/mnt/weka/mbedrosian/code/casf-benchmark')
W=Path('/mnt/weka/mbedrosian/pharma_generation_analysis/v2_1p7b_coordtoken_step41054')
OLD='qwen_1p7b_fsq_bigdata_step47023';NEW='qwen_1p7b_v2_coordtoken_step41054';NAME='Qwen 1.7B FSQ'
DB=B/'data/results/casf_analysis_dashboard.sqlite'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def edit_csv(path,fn):
 with path.open() as f:r=csv.DictReader(f);fields=r.fieldnames;rows=list(r)
 for row in rows:fn(row)
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
def replace_records(path,source,col,identity,fields=None):
 # Keep all non-Qwen numerical strings, including their precision, unchanged.
 source=source.set_index(identity);source.index=source.index.map(str)
 def change(r):
  if not r[col].startswith((OLD,NEW)):return
  key=r[identity]
  if identity=='method':key=key.replace(OLD,NEW)
  s=source.loc[key];assert isinstance(s,pd.Series)
  for k in (fields or r.keys()):
   if k==identity:r[k]=key
   elif k in s.index:r[k]='' if pd.isna(s[k]) else str(s[k])
 edit_csv(path,change)
with sqlite3.connect(f'file:{DB}?mode=ro',uri=True) as c:
 casf=pd.read_sql_query('SELECT * FROM per_ligand_long WHERE family=?',c,params=[NEW])
 base=pd.read_sql_query("SELECT * FROM per_ligand_long WHERE method='chembl3d_gt_pb' AND ligand_set IN ('core','ref_dedup')",c)
casf.to_csv(HERE/'casf-per-entry.csv',index=False)
PL=W/'plinder_panel'/NEW
pl=pd.read_csv(PL/'per_molecule.csv');summary=pd.read_csv(PL/'summary.csv')
assert len(pl)==pl.smiles.nunique()==22 and pl.num_true_confs.sum()==1067
for filename in ['preflight.json','submission.json']:
 shutil.copy2(W/filename,HERE/filename)
shutil.copy2(W/'gen_results'/f'{NEW}_plinder'/'generation_metadata.json',HERE/'plinder-generation-metadata.json')
for number in [1,2,3]:
 def change(r):
  if r['method']!=NAME:return
  cohort=r['cohort']
  if cohort=='plinder':
   d=pl;distance=d.mat_r;clusters=d.greedy_clusters_0p5
   hit5=d.cov_r_0p5.mean()*100;hit75=d.cov_r_0p75.mean()*100
  else:
   tier=r['budget'] if number==3 else 'chembl_count'
   d=casf[casf.ligand_set.eq('core' if cohort=='core' else 'ref_dedup')&casf.tier.eq(tier)]
   distance=d.casf_best_rmsd;clusters=d.greedy_clusters_0p5
   hit5=100*distance.le(.5).mean();hit75=100*distance.le(.75).mean()
  if number==1:r.update(n_defined=str(clusters.count()),clusters_0p5=str(clusters.mean()),clusters_1p0=str(d.greedy_clusters_1p0.mean()),largest_percent=str(100*d.largest_cluster_fraction_1p0.mean()),status=NEW)
  elif number==2:r.update(n_distance=str(distance.count()),recovery_0p5=str(hit5),recovery_0p75=str(hit75),distance=str(distance.mean()))
  else:r.update(n_defined=str(distance.count()),distance=str(distance.mean()),recovery=str(hit75),clusters=str(clusters.mean()),mean_retained=str(d.post_pb_confs.fillna(0).mean()),status=NEW)
 edit_csv(P/f'results-data/three-dataset-tables/table{number}.csv',change)
# Current single-reference entry archives.
for path,cohort in [(P/'results-data/ref-panel/per-entry.csv','ref_dedup'),(P/'figure-data/candidates/core-entries.csv','core')]:
 lookup=casf[casf.ligand_set.eq(cohort)].set_index(['method','mol_id'])
 def change(r):
  if not r['method'].startswith((OLD,NEW)):return
  r['method']=r['method'].replace(OLD,NEW);s=lookup.loc[(r['method'],r['mol_id'])]
  for k in r:
   if k in s.index:r[k]='' if pd.isna(s[k]) else str(s[k])
 edit_csv(path,change)
# Current PLINDER summaries, per-molecule values, and all nearest distances.
for path in [P/'plinder22-data/per-molecule.csv',P/'results-data/three-dataset-tables/plinder-per-molecule.csv']:
 replace_records(path,pl,'label','smiles')
replace_records(P/'results-data/druglike_selected.csv',pl.rename(columns={'cov_r_0p75':'cov_r_075','cov_p_0p75':'cov_p_075'}),'label','smiles')
def sumchange(r):
 if r['label'] in [OLD,NEW]:
  for k in r:
   if k in summary:r[k]=str(summary.iloc[0][k])
edit_csv(P/'plinder22-data/summary.csv',sumchange)
entries=json.loads((P/'plinder22-data/entries.json').read_text());oldz=np.load(P/'plinder22-data/nearest-distances.npz');minima={k:oldz[k] for k in oldz.files}
z=np.load(W/'plinder_cache'/f'{NEW}.npz',allow_pickle=True);mats={str(s):z[f'm{i}'] for i,s in enumerate(z['smiles'])}
curves=[];ts=np.linspace(0,1.5,301);by=pl.set_index('smiles');recalls=[];precisions=[]
for e in entries:
 if e['method'] not in [OLD,NEW]:continue
 e['method']=NEW;m=np.asarray(mats[e['smiles']],dtype=float);r=m.min(axis=1);p=m.min(axis=0);row=by.loc[e['smiles']]
 assert m.shape==(row.num_true_confs,row.n_generated) and np.isfinite(m).all()
 np.testing.assert_allclose([r.mean(),p.mean()],[row.mat_r,row.mat_p],atol=1e-12,rtol=0)
 for t,tag in [(.5,'0p5'),(.75,'0p75'),(1.,'1p0')]:np.testing.assert_allclose([(r<=t).mean(),(p<=t).mean()],[row['cov_r_'+tag],row['cov_p_'+tag]],atol=1e-12,rtol=0)
 minima[e['key']+'_reference']=r;minima[e['key']+'_generated']=p
 recalls.append((r[:,None]<=ts).mean(axis=0));precisions.append((p[:,None]<=ts).mean(axis=0))
np.savez_compressed(P/'plinder22-data/nearest-distances.npz',**minima)
(P/'plinder22-data/entries.json').write_text(json.dumps(entries,indent=2)+'\n')
rc=np.mean(recalls,axis=0);pc=np.mean(precisions,axis=0)
def curvechange(r):
 if r['method'] in [OLD,NEW]:
  i=round(float(r['threshold'])/.005);r.update(method=NEW,recall=str(rc[i]),precision=str(pc[i]))
edit_csv(P/'plinder22-data/threshold-curves.csv',curvechange)
# Common baseline descriptors preserve the existing difficulty bins.
for filename,descriptor,labels,bins in [('casf-size-source.csv','heavy_atoms',['<20','20-29','30-39','40+'],[0,20,30,40,np.inf]),('casf-flexibility-source.csv','rotors',['0-3','4-6','7-8','9+'],[-1,3,6,8,np.inf])]:
 def change(r):
  if not r['method'].startswith((OLD,NEW)):return
  c='core' if r['set']=='core' else 'ref_dedup';b=base[base.ligand_set.eq(c)].set_index('mol_id')
  g=pd.cut(b.heavy_atoms if descriptor=='heavy_atoms' else b.rotatable_bonds,bins,labels=labels,right=descriptor!='heavy_atoms')
  d=casf[casf.ligand_set.eq(c)&casf.tier.eq('fixed')].set_index('mol_id').loc[b.index];d=d[g.eq(r[descriptor])]
  hits=int(d.casf_best_rmsd.le(.75).sum());r.update(method=NEW+'_fixed',n=str(len(d)),hit_075=str(hits/len(d)) if len(d) else '')
  if 'hits' in r:r['hits']=str(hits)
  if 'rmsd' in r:r['rmsd']=str(d.casf_best_rmsd.mean())
 edit_csv(P/'results-data/size-flexibility'/filename,change)
strata={}
for f,desc,col in [('casf-size-source.csv','size','heavy_atoms'),('casf-flexibility-source.csv','rotors','rotors')]:
 d=pd.read_csv(P/'results-data/size-flexibility'/f)
 for r in d[d['set'].eq('core')&d.method.eq(NEW+'_fixed')].to_dict('records'):strata[desc,r[col]]=r
def stratac(r):
 if r['key'] not in [OLD,NEW]:return
 s=strata[r['descriptor'],r['group']];r.update(key=NEW,n=str(s['n']),hits=str(round(s['n']*s['hit_075'])),hit=str(100*s['hit_075']))
edit_csv(P/'results-data/strata.csv',stratac)
def points(r):
 if not r['method'].startswith((OLD,NEW)):return
 d=casf[casf.ligand_set.eq('core')&casf.tier.eq('chembl_count')]
 r.update(method=NEW+'_chembl_count',clusters_0p5=str(d.greedy_clusters_0p5.mean()),best_rmsd=str(d.casf_best_rmsd.mean()),n_defined=str(d.casf_best_rmsd.count()))
edit_csv(P/'results-data/diversity-rmsd-points.csv',points)
def tierchange(r):
 if r['key'] not in [OLD,NEW] or r['cohort']!='core':return
 d=casf[casf.ligand_set.eq('core')&casf.tier.eq(r['tier'])]
 r.update(key=NEW,evaluated=str(d.casf_best_rmsd.count()),mean_retained=str(d.post_pb_confs.fillna(0).mean()),hit_0p75=str(100*d.casf_best_rmsd.le(.75).mean()),mean_minimum_rmsd=str(d.casf_best_rmsd.mean()),mean_ensemble_median_rmsd=str(d.casf_median_rmsd.mean()))
edit_csv(P/'figure-data/recovery-table-tiers.csv',tierchange)
# Figure 6 keeps the separate 0.6B run and its nine-method scope.
track=P/'plinder22-data/molecule-tracking'
replace_records(track/'molecule_matching_values.csv',pl,'label','smiles')
methods=json.loads((track/'methods.json').read_text());methods={k.replace(OLD,NEW):v for k,v in methods.items()}
d=pd.read_csv(track/'molecule_matching_values.csv');order=d.groupby('label').mat_r.mean().sort_values().index
(track/'methods.json').write_text(json.dumps({k:methods[k] for k in order},indent=2)+'\n')
means=d.groupby('label',sort=False)[['mat_r','mat_p']].mean().loc[order];means.to_csv(track/'matching_means.csv')
# Exact identities in active renderers/configuration. Historical archives stay intact.
paths=['sampling-data/palette.json','sampling-data/display.json','plinder22-data/methods.json','plinder22-data/plot.py','results-data/build_druglike_table.py','build_figures.py','results-data/ref-panel/build_supporting.py','plinder22-data/molecule-tracking/plot.py']
for rel in paths:
 p=P/rel;p.write_text(p.read_text().replace(OLD,NEW))
prov=P/'figure-data/candidates/provenance.json';v=json.loads(prov.read_text());v['export_sha256']=sha(P/'figure-data/candidates/core-entries.csv')
for s in v['sources']:s['sha256']=sha(P/'figure-data/candidates'/s['path'])
v['primary_qwen_update']='results-data/qwen-v2/provenance.json';prov.write_text(json.dumps(v,indent=2)+'\n')
inputs=[DB,PL/'per_molecule.csv',PL/'summary.csv',W/'plinder_cache'/f'{NEW}.npz',W/'preflight.json']
(HERE/'provenance.json').write_text(json.dumps({'checkpoint':NEW,'replaces':OLD,'display_label':NAME,'date':'2026-11-04','cohorts':{'core':94,'ref':1044,'plinder':22},'plinder_references':1067,'source_sha256':{str(p):sha(p) for p in inputs},'scope':'Replace primary Qwen only; other methods and the separate Qwen 0.6B CASF-dedup result are preserved. Energy and sampling protocol measurements are separately archived in calculation-provenance.json.','training_exclusions':'Not inferred from the checkpoint label; executed corpus/exclusion records still require documentation.'},indent=2)+'\n')
print('Updated primary Qwen table inputs, PLINDER distances, strata and review-figure inputs.')
