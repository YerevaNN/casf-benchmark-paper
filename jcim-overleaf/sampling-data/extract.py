"""Extract aligned RMSDs from the exact retained fixed-target conformers.
Run with PYTHONPATH=src in the analysis checkout. No geometry is optimized.
"""
from pathlib import Path
import concurrent.futures as futures
import hashlib,json,sqlite3
import numpy as np
import pandas as pd
from rdkit import Chem,RDLogger,rdBase
from casf_benchmark.analysis.metrics import best_aligned_rmsd
from casf_benchmark.cli.analyze_conformer_sets import load_casf_ligand
from casf_benchmark.paths import DEFAULT_CORE_LIGAND_DIR
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
CONFIG=json.loads((OUT/'config.json').read_text());CACHE=OUT/'cache'
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def task(job):
    RDLogger.DisableLog('rdApp.*')
    path=Path(job['sdf']);digest=sha(path) if path.exists() else None
    fingerprint=hashlib.sha256(json.dumps(dict(job=job,sdf_sha256=digest,script=sha(Path(__file__)),rdkit=rdBase.rdkitVersion),sort_keys=True).encode()).hexdigest()
    cache=CACHE/(job['method']+'__'+job['mol_id']+'.json')
    if cache.exists():
        old=json.loads(cache.read_text())
        if old['fingerprint']==fingerprint:return old
    mols=[m for m in Chem.SDMolSupplier(str(path),removeHs=False)] if path.exists() and path.stat().st_size else []
    assert len(mols)==job['n_retained'] and all(m is not None for m in mols),(job['method'],job['mol_id'],len(mols),job['n_retained'])
    ref=load_casf_ligand(job['mol_id'],DEFAULT_CORE_LIGAND_DIR)
    assert ref is not None
    values=np.array([best_aligned_rmsd(m,ref) for m in mols])
    finite=values[np.isfinite(values)]
    if job['best'] is not None:
        assert len(finite)>0 and np.isclose(finite.min(),job['best'],rtol=1e-5,atol=1e-5),(job['method'],job['mol_id'],finite.min(),job['best'])
        assert np.isclose(np.median(finite),job['median'],rtol=1e-5,atol=1e-5)
    else:assert len(finite)==0
    result=dict(fingerprint=fingerprint,method=job['method'],mol_id=job['mol_id'],n_candidates=job['n_candidates'],n_retained=len(mols),sdf=str(path),sdf_sha256=digest,
        rmsd=[float(x) if np.isfinite(x) else None for x in values],n_rmsd_failed=int((~np.isfinite(values)).sum()))
    cache.write_text(json.dumps(result)+'\n');return result

def main():
    CACHE.mkdir(exist_ok=True)
    db=ROOT/'data/results/casf_analysis_dashboard.sqlite';assert sha(db)==CONFIG['database_sha256']
    with sqlite3.connect(f'file:{db}?mode=ro',uri=True) as c:
        per=pd.read_sql_query('SELECT * FROM per_ligand_long WHERE ligand_set="core"',c)
        sources=pd.read_sql_query('SELECT * FROM analysis_sources',c).set_index('run_id').root.to_dict()
    selected=pd.read_csv(ROOT/'docs/results_tables/casf_selected.csv')
    methods=selected[selected.tier.eq('fixed')].method.tolist();ids=sorted(per[per.method.eq('chembl3d_gt_pb')].mol_id)
    assert len(ids)==94
    jobs=[]
    for method in methods:
        for r in per[per.method.eq(method)].sort_values('mol_id').to_dict('records'):
            assert 0<=r['post_pb_confs']<=r['pb_input_confs']<=1000
            jobs.append(dict(method=method,mol_id=r['mol_id'],sdf=str(Path(sources[r['run_id']])/'generation'/method/(r['mol_id']+'.sdf')),
                n_candidates=int(r['pb_input_confs']),n_retained=int(r['post_pb_confs']),best=float(r['casf_best_rmsd']) if pd.notna(r['casf_best_rmsd']) else None,median=float(r['casf_median_rmsd']) if pd.notna(r['casf_median_rmsd']) else None))
    records=[]
    with futures.ProcessPoolExecutor(max_workers=CONFIG['workers']) as pool:
        for i,r in enumerate(pool.map(task,jobs),1):
            records.append(r)
            if i%50==0:print(f'{i}/{len(jobs)} pools checked',flush=True)
    # Rejected candidates and RMSD failures are non-hits. Their ordering is irrelevant
    # for uniform subsets. Padding beyond N is not part of the sampling population.
    distances=np.full((len(methods),94,1000),np.inf);n_candidates=np.zeros((len(methods),94),int);audit=[]
    for r in records:
        m=methods.index(r['method']);i=ids.index(r['mol_id']);n_candidates[m,i]=r['n_candidates']
        distances[m,i,:r['n_retained']]=[x if x is not None else np.inf for x in r['rmsd']]
        audit.append({k:v for k,v in r.items() if k!='rmsd'})
    distances.sort(axis=2)
    ref=per[per.method.eq('chembl3d_gt_pb')].set_index('mol_id').reindex(ids).casf_best_rmsd.to_numpy()
    np.savez_compressed(OUT/'candidate_rmsds.npz',rmsd=distances,n_candidates=n_candidates,methods=np.array(methods),mol_ids=np.array(ids),reference_best_rmsd=ref)
    pd.DataFrame(audit).to_csv(OUT/'pool_audit.csv',index=False)
    provenance=dict(database=str(db),database_sha256=sha(db),rdkit=rdBase.rdkitVersion,extract_sha256=sha(Path(__file__)),config_sha256=sha(OUT/'config.json'),metric_sha256=sha(ROOT/'src/casf_benchmark/analysis/metrics.py'),ligand_directory=str(DEFAULT_CORE_LIGAND_DIR),n_conformers=sum(r['n_retained'] for r in records),n_rmsd_failed=sum(r['n_rmsd_failed'] for r in records),protocol='Uniform subsets of archived pre-PoseBusters candidate pools; k capped at available candidates; rejected candidates and unmeasurable RMSDs count as non-hits. Budget is not raw generator attempts or computing cost.')
    (OUT/'extraction-provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(provenance,indent=2),flush=True)
if __name__=='__main__':main()
