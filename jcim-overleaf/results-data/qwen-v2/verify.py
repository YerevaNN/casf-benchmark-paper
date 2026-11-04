"""Check the primary-checkpoint substitution independently of the renderers."""
from pathlib import Path
import hashlib,io,json,re,sqlite3,subprocess
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
H=Path(__file__).resolve().parent;P=H.parents[1];REPO=P.parent
NEW='qwen_1p7b_v2_coordtoken_step41054';OLD='qwen_1p7b_fsq_bigdata_step47023'
prov=json.loads((H/'provenance.json').read_text());database=next(p for p in prov['source_sha256'] if p.endswith('casf_analysis_dashboard.sqlite'))
with sqlite3.connect(f'file:{database}?mode=ro',uri=True) as c:
 source=pd.read_sql_query('SELECT * FROM per_ligand_long WHERE family=?',c,params=[NEW])
arch=pd.read_csv(H/'casf-per-entry.csv');pd.testing.assert_frame_equal(source,arch,check_dtype=False,check_exact=False)
for cohort,n in [('core',94),('ref',1044)]:
 d=source[source.ligand_set.eq('core' if cohort=='core' else 'ref_dedup')]
 assert d.groupby('tier').size().eq(n).all()
 for number in [1,2,3]:
  table=pd.read_csv(P/f'results-data/three-dataset-tables/table{number}.csv')
  for r in table[table.cohort.eq(cohort)&table.method.eq('Qwen 1.7B FSQ')].itertuples():
   q=d[d.tier.eq(r.budget if number==3 else 'chembl_count')]
   if number==1:a=[r.n_defined,r.clusters_0p5,r.clusters_1p0,r.largest_percent];b=[q.greedy_clusters_0p5.count(),q.greedy_clusters_0p5.mean(),q.greedy_clusters_1p0.mean(),100*q.largest_cluster_fraction_1p0.mean()]
   elif number==2:a=[r.n_distance,r.distance,r.recovery_0p5,r.recovery_0p75];b=[q.casf_best_rmsd.count(),q.casf_best_rmsd.mean(),100*q.casf_best_rmsd.le(.5).mean(),100*q.casf_best_rmsd.le(.75).mean()]
   else:a=[r.n_defined,r.distance,r.recovery,r.clusters,r.mean_retained];b=[q.casf_best_rmsd.count(),q.casf_best_rmsd.mean(),100*q.casf_best_rmsd.le(.75).mean(),q.greedy_clusters_0p5.mean(),q.post_pb_confs.fillna(0).mean()]
   np.testing.assert_allclose(a,b,rtol=0,atol=1e-10)
 z=np.load(H/cohort/'candidate_rmsds.npz');dist=z['rmsd'][0];counts=z['n_candidates'][0];ids=z['mol_ids']
 q=d[d.tier.eq('fixed')].set_index('mol_id').loc[ids]
 np.testing.assert_array_equal(counts,q.pb_input_confs)
 np.testing.assert_array_equal(np.isfinite(dist).sum(axis=1),q.post_pb_confs)
 np.testing.assert_allclose(dist.min(axis=1),q.casf_best_rmsd,atol=1e-5,rtol=1e-5)
 curves=pd.read_csv(H/cohort/'recovery_curves.csv')
 for r in curves.itertuples():
  expected=hypergeom.sf(0,counts,(dist<=r.threshold).sum(axis=1),np.minimum(counts,r.budget)).mean()
  np.testing.assert_allclose(r.recovery,expected,atol=1e-12,rtol=0)
 e=pd.read_csv(H/cohort/'energy-per-conformer.csv.gz');per=pd.read_csv(H/cohort/'energy-per-molecule.csv').set_index('mol_id')
 assert e.max_heavy_coordinate_change.max()<=1e-10 and e.status.eq('evaluated').all()
 stats=e.groupby('mol_id').energy.agg(['count','mean','min','max','median'])
 for col,src in [('count','n_energy'),('mean','mean_energy'),('min','min_energy'),('max','max_energy'),('median','median_energy')]:np.testing.assert_allclose(stats[col],per.loc[stats.index,src],rtol=1e-12,atol=1e-9)
 np.testing.assert_allclose(e.groupby('mol_id').energy.std(ddof=0),per.loc[stats.index,'std_energy'],rtol=1e-12,atol=1e-9)
# Other observations remain identical except the subsequently approved independent
# PLINDER torsion replacement, checked by plinder22-data/independent-torsion/verify.py.
for rel,key in [('plinder22-data/per-molecule.csv','label'),('results-data/ref-panel/per-entry.csv','method'),('figure-data/candidates/core-entries.csv','method'),('energy-data/per_molecule.csv','method'),('results-data/core-ref-protocol/ref/per_molecule.csv','method'),('sampling-data/recovery_curves.csv','method'),('results-data/core-ref-protocol/ref/recovery_curves.csv','method')]:
 old=pd.read_csv(io.StringIO(subprocess.check_output(['git','show',f'ed8b4bf:jcim-overleaf/{rel}'],cwd=REPO,text=True)))
 new=pd.read_csv(P/rel)
 old=old[~old[key].str.startswith(OLD)].reset_index(drop=True);new=new[~new[key].str.startswith(NEW)].reset_index(drop=True)
 if rel=='plinder22-data/per-molecule.csv':
  old=old[~old[key].isin(['torsion_raw_druglike','torsion_minimized_22_druglike'])].reset_index(drop=True)
  new=new[~new[key].isin(['torsion_independent_raw_druglike','torsion_independent_minimized_druglike'])].reset_index(drop=True)
 pd.testing.assert_frame_equal(old,new,check_exact=True)
# Check highlighted ranks from unrounded values, excluding the stored baseline.
for number,name,columns in [(1,'diversity',['clusters_0p5','clusters_1p0','largest_percent']),(2,'count-recovery',['recovery_0p5','recovery_0p75','distance'])]:
 df=pd.read_csv(P/f'results-data/three-dataset-tables/table{number}.csv');tex=(P/f'tables/{name}.tex').read_text()
 for cohort_index,cohort in enumerate(['plinder','core','ref']):
  peers=df[df.cohort.eq(cohort)&df.method.ne('ChEMBL3D-PB')]
  for metric_index,metric in enumerate(columns):
   ranks=sorted(set(peers[metric].dropna()),reverse=metric not in ['distance','largest_percent'])[:2]
   for r in peers.to_dict('records'):
    line=next(l for l in tex.splitlines() if l.startswith(r['method']+' &') or l.startswith(r['method']+'$'))
    cell=line.split(' & ')[1+cohort_index*3+metric_index]
    assert (r'\textbf{' in cell)==(r[metric]==ranks[0]),(name,cohort,r['method'],metric)
    assert (r'\underline{' in cell)==(len(ranks)>1 and r[metric]==ranks[1]),(name,cohort,r['method'],metric)
# Old primary checkpoint cannot remain in active model configurations or manuscript prose.
for rel in ['sampling-data/palette.json','sampling-data/display.json','plinder22-data/methods.json','plinder22-data/molecule-tracking/methods.json','results.tex','methods.tex','tables/generation-methods.tex']:
 assert OLD not in (P/rel).read_text() and '47,023' not in (P/rel).read_text(),rel
print('Verified exact checkpoint identity, all CASF replacement aggregates, finite-pool probabilities, unoptimized energy summaries, unchanged comparator observations, and best/second-best table highlighting.')
