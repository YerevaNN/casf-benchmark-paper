"""Measure the replacement checkpoint under the existing paper protocols.

Saved PB-passing coordinates only; no optimization or validity replay.
Run with the benchmark src on PYTHONPATH. Cache is outside the paper repository.
"""
from pathlib import Path
import concurrent.futures as futures
import hashlib, importlib.util, json, sqlite3, sys
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger, rdBase
from casf_benchmark.analysis.metrics import best_aligned_rmsd
from casf_benchmark.cli.analyze_conformer_sets import load_casf_ligand
from casf_benchmark.paths import DEFAULT_CORE_LIGAND_DIR, DEFAULT_CASF_REF_LIGAND_DIR
HERE=Path(__file__).resolve().parent
PAPER=HERE.parents[1]
BENCH=Path('/mnt/weka/mbedrosian/code/casf-benchmark')
WORK=Path('/mnt/weka/mbedrosian/pharma_generation_analysis/v2_1p7b_coordtoken_step41054')
NEW='qwen_1p7b_v2_coordtoken_step41054'
CACHE=WORK/'paper_protocol_cache';CACHE.mkdir(exist_ok=True)
sys.path.insert(0,str(PAPER/'energy-data'))
from rescore import rescore

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def task(r):
 RDLogger.DisableLog('rdApp.*')
 cohort='core' if r['ligand_set']=='core' else 'ref'
 path=WORK/'casf_benchmark_runs'/f'{NEW}_{cohort}'/'generation'/r['method']/(r['mol_id']+'.sdf')
 fingerprint={'sdf_sha256':sha(path),'protocol':'paper-common-H-MMFF94s-and-symmetry-RMSD-v1','rdkit':rdBase.rdkitVersion}
 cache=CACHE/f"{cohort}_{r['tier']}_{r['mol_id']}.json"
 if cache.exists():
  d=json.loads(cache.read_text())
  if d['fingerprint']==fingerprint:return d
 mols=[m for m in Chem.SDMolSupplier(str(path),removeHs=False)] if path.stat().st_size else []
 assert len(mols)==r['post_pb_confs'] and all(m is not None for m in mols),(str(path),len(mols),r['post_pb_confs'],sum(m is None for m in mols))
 if r['tier']=='fixed':
  ref=load_casf_ligand(r['mol_id'],DEFAULT_CORE_LIGAND_DIR if cohort=='core' else DEFAULT_CASF_REF_LIGAND_DIR)
  assert ref is not None
  vals=np.array([best_aligned_rmsd(m,ref) for m in mols]);finite=vals[np.isfinite(vals)]
  assert np.isclose(finite.min(),r['casf_best_rmsd'],atol=1e-5,rtol=1e-5),(r['mol_id'],finite.min(),r['casf_best_rmsd'])
  assert np.isclose(np.median(finite),r['casf_median_rmsd'],atol=1e-5,rtol=1e-5)
  values=[float(v) if np.isfinite(v) else None for v in vals]
 else:
  values=[rescore(m) for m in mols]
 d={'fingerprint':fingerprint,'sdf':str(path),'cohort':cohort,'mol_id':r['mol_id'],'tier':r['tier'],'values':values}
 cache.write_text(json.dumps(d)+'\n');return d

def main():
 db=BENCH/'data/results/casf_analysis_dashboard.sqlite'
 with sqlite3.connect(f'file:{db}?mode=ro',uri=True) as c:
  rows=pd.read_sql_query('SELECT * FROM per_ligand_long WHERE family=?',c,params=[NEW])
  bases=pd.read_sql_query("SELECT * FROM per_ligand_long WHERE method='chembl3d_gt_pb' AND ligand_set IN ('core','ref_dedup')",c)
 rows.to_csv(HERE/'casf-per-entry.csv',index=False)
 jobs=rows[rows.tier.isin(['fixed','chembl_count'])].to_dict('records')
 results=[]
 with futures.ProcessPoolExecutor(max_workers=8) as pool:
  for i,d in enumerate(pool.map(task,jobs),1):
   results.append(d)
   if i%50==0:print(f'{i}/{len(jobs)} saved ensembles measured',flush=True)
 for cohort,n in [('core',94),('ref',1044)]:
  dest=HERE/cohort;dest.mkdir(exist_ok=True)
  sub=rows[rows.ligand_set.eq('core' if cohort=='core' else 'ref_dedup')]
  ids=sorted(sub.mol_id.unique());assert len(ids)==n
  by={(r['tier'],r['mol_id']):r for r in sub.to_dict('records')}
  energy=[];records=[];dist=np.full((1,n,1000),np.inf);counts=np.zeros((1,n),int)
  for d in results:
   if d['cohort']!=cohort:continue
   r=by[d['tier'],d['mol_id']]
   if d['tier']=='fixed':
    i=ids.index(d['mol_id']);counts[0,i]=r['pb_input_confs']
    values=[v if v is not None else np.inf for v in d['values']]
    dist[0,i,:len(values)]=values
   else:
    vals=np.array([v['energy'] for v in d['values'] if v['energy'] is not None]);assert len(vals)==len(d['values'])
    energy.append(dict(method=NEW+'_chembl_count',mol_id=d['mol_id'],n_retained=len(d['values']),n_energy=len(vals),n_failed=len(d['values'])-len(vals),mean_energy=vals.mean() if len(vals) else np.nan,std_energy=vals.std() if len(vals) else np.nan,min_energy=vals.min() if len(vals) else np.nan,max_energy=vals.max() if len(vals) else np.nan,median_energy=np.median(vals) if len(vals) else np.nan,clusters_0p5=r['greedy_clusters_0p5'],archived_std=r['energy_std']))
    records.extend(dict(method=NEW+'_chembl_count',mol_id=d['mol_id'],conformer_index=i,**v) for i,v in enumerate(d['values']))
  pd.DataFrame(energy).to_csv(dest/'energy-per-molecule.csv',index=False)
  pd.DataFrame(records).to_csv(dest/'energy-per-conformer.csv.gz',index=False)
  dist.sort(axis=2)
  refs=bases[bases.ligand_set.eq('core' if cohort=='core' else 'ref_dedup')].set_index('mol_id').loc[ids].casf_best_rmsd.to_numpy()
  np.savez_compressed(dest/'candidate_rmsds.npz',rmsd=dist,n_candidates=counts,methods=np.array([NEW+'_fixed']),mol_ids=np.array(ids),reference_best_rmsd=refs)
  spec=importlib.util.spec_from_file_location('sampling',PAPER/'sampling-data/analyze.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  module.OUT=dest;module.main()
 audit=[{k:v for k,v in d.items() if k!='values'} for d in results]
 (HERE/'calculation-provenance.json').write_text(json.dumps({'checkpoint':NEW,'database':str(db),'database_sha256':sha(db),'rdkit':rdBase.rdkitVersion,'energy_protocol':json.loads((PAPER/'energy-data/protocol.json').read_text()),'sampling_protocol':json.loads((PAPER/'sampling-data/config.json').read_text()),'inputs':audit,'coordinate_optimization':False},indent=2)+'\n')
if __name__=='__main__':main()
