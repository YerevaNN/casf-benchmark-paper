"""Rescore retained ChEMBL-count ensembles with a common hydrogen protocol.
Run in the analysis checkout: python docs/energy_analysis/rescore.py
No heavy-atom coordinates, source molecules, or benchmark databases are written.
"""
from pathlib import Path
import concurrent.futures as futures
import csv
import hashlib
import json
import math
import os
import sqlite3
import gzip
import numpy as np
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
CONFIG = json.loads((OUT/'protocol.json').read_text())
CACHE = OUT/'cache'
DB = ROOT/'data/results/casf_analysis_dashboard.sqlite'

def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def rescore(mol):
    heavy = Chem.RemoveHs(mol)
    xyz = np.asarray(heavy.GetConformer().GetPositions()).copy()
    if not np.isfinite(xyz).all():
        return dict(energy=None, status='nonfinite_coordinates', original_hydrogens=0)
    original_h = sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
    prepared = Chem.AddHs(heavy, addCoords=True)
    props = AllChem.MMFFGetMoleculeProperties(prepared, mmffVariant=CONFIG['force_field'])
    if props is None:
        return dict(energy=None, status='missing_parameters', original_hydrogens=original_h)
    ff = AllChem.MMFFGetMoleculeForceField(prepared, props,
        nonBondedThresh=CONFIG['nonbonded_threshold'],
        ignoreInterfragInteractions=CONFIG['ignore_interfragment_interactions'])
    for i in range(heavy.GetNumAtoms()):
        ff.AddFixedPoint(i)
    ff.Initialize()
    status = 1
    for _ in range(CONFIG['max_minimization_calls']):
        status = ff.Minimize(maxIts=CONFIG['max_iterations'],
            forceTol=CONFIG['force_tolerance'], energyTol=CONFIG['energy_tolerance'])
        if status == 0:
            break
    change = np.max(np.abs(np.asarray(prepared.GetConformer().GetPositions())[:len(xyz)]-xyz))
    assert change <= CONFIG['heavy_coordinate_tolerance_angstrom'], change
    energy = float(ff.CalcEnergy())
    return dict(energy=energy if math.isfinite(energy) else None,
        status='converged' if status == 0 and math.isfinite(energy) else 'not_converged',
        original_hydrogens=original_h, explicit_hydrogens=prepared.GetNumAtoms()-heavy.GetNumAtoms(),
        max_heavy_coordinate_change=float(change))

def task(job):
    RDLogger.DisableLog('rdApp.*')
    cache = CACHE/(job['method']+'__'+job['mol_id']+'.json')
    job_hash = hashlib.sha256(json.dumps(job,sort_keys=True).encode()).hexdigest()
    if cache.exists():
        previous=json.loads(cache.read_text())
        if previous['job_sha256']==job_hash:
            return previous
    if job['method']=='chembl3d_gt_pb':
        import pandas as pd
        from casf_benchmark.cli.analyze_conformer_sets import get_reference_chembl_mols
        from casf_benchmark.paths import DEFAULT_CHEMBL_DATASET_ROOT
        mols,_=get_reference_chembl_mols(pd.Series(job['mapping']),job['mol_id'],
            DEFAULT_CHEMBL_DATASET_ROOT/'topologies',DEFAULT_CHEMBL_DATASET_ROOT/'zarr_database')
        # This core reference snapshot has zero PB rejections; verify count and energy fingerprint.
        from casf_benchmark.analysis.metrics import energy_stats
        old=energy_stats(mols)
        for key in ['energy_median','energy_std']:
            assert math.isclose(old[key],job[key],rel_tol=1e-5,abs_tol=1e-4),(job['mol_id'],key,old[key],job[key])
    else:
        path=Path(job['sdf_path'])
        mols=list(Chem.SDMolSupplier(str(path),removeHs=False)) if path.exists() and path.stat().st_size else []
        assert all(m is not None for m in mols),job['sdf_path']
    assert len(mols)==job['post_pb_confs'],(job['method'],job['mol_id'],len(mols),job['post_pb_confs'])
    records=[]
    for i,mol in enumerate(mols):
        try:
            result=rescore(mol)
        except AssertionError:
            raise
        except Exception as exc:
            result=dict(energy=None,status=type(exc).__name__+': '+str(exc))
        records.append(dict(conformer_index=i,**result))
    result=dict(job_sha256=job_hash,method=job['method'],mol_id=job['mol_id'],records=records)
    cache.write_text(json.dumps(result)+'\n')
    return result

def main():
    assert digest(DB)==CONFIG['source_database_sha256']
    CACHE.mkdir(exist_ok=True)
    con=sqlite3.connect('file:'+str(DB)+'?mode=ro',uri=True);con.row_factory=sqlite3.Row
    sources={x['run_id']:Path(x['root']) for x in con.execute('select * from analysis_sources')}
    labels=list(csv.DictReader((ROOT/'docs/results_tables/clustering_radius_comparison.csv').open()))
    from casf_benchmark.paths import DEFAULT_CHEMBL_MAP_CSV
    mapping={r['ligand_id']:r for r in csv.DictReader(DEFAULT_CHEMBL_MAP_CSV.open())}
    fingerprint=dict(config_sha256=digest(OUT/'protocol.json'),script_sha256=digest(__file__),
        mapping_sha256=digest(DEFAULT_CHEMBL_MAP_CSV),rdkit=rdBase.rdkitVersion)
    jobs=[];archive=[];input_files={}
    for method in labels:
        for r in con.execute('select * from per_ligand_long where ligand_set=? and method=?',('core',method['method'])):
            j=dict(method=r['method'],mol_id=r['mol_id'],post_pb_confs=int(r['post_pb_confs']),fingerprint=fingerprint)
            archive.append(dict(r))
            if r['method']=='chembl3d_gt_pb':
                raw=con.execute('select pb_fail_confs from per_ligand_long where ligand_set=? and method=? and mol_id=?',('core','chembl3d_gt',r['mol_id'])).fetchone()
                assert raw[0]==0
                j.update(mapping=mapping[r['mol_id']],energy_median=r['energy_median'],energy_std=r['energy_std'])
            else:
                path=sources[r['run_id']]/'generation'/r['method']/(r['mol_id']+'.sdf')
                j['sdf_path']=str(path);j['sdf_sha256']=digest(path) if path.exists() else None
                if r['post_pb_confs']:assert path.exists(),path
                input_files[str(path)]=j['sdf_sha256']
            jobs.append(j)
    results=[]
    with futures.ProcessPoolExecutor(max_workers=CONFIG['workers']) as pool:
        for i,result in enumerate(pool.map(task,jobs),1):
            results.append(result)
            if i%50==0:print(f'{i}/{len(jobs)} ensembles rescored',flush=True)
    with gzip.open(OUT/'per_conformer.csv.gz','wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['method','mol_id','conformer_index','energy','status','original_hydrogens','explicit_hydrogens','max_heavy_coordinate_change']);w.writeheader()
        for result in results:
            for record in result['records']:w.writerow(dict(method=result['method'],mol_id=result['mol_id'],**record))
    by_id={(r['method'],r['mol_id']):r for r in archive}
    per=[]
    for result in results:
        old=by_id[result['method'],result['mol_id']]
        vals=np.array([x['energy'] for x in result['records'] if x['energy'] is not None],float)
        per.append(dict(method=result['method'],mol_id=result['mol_id'],n_retained=old['post_pb_confs'],n_energy=len(vals),
            n_unconverged=sum(x['status']!='converged' for x in result['records']),
            mean_energy=float(vals.mean()) if len(vals) else None,std_energy=float(vals.std(ddof=0)) if len(vals) else None,
            min_energy=float(vals.min()) if len(vals) else None,max_energy=float(vals.max()) if len(vals) else None,
            median_energy=float(np.median(vals)) if len(vals) else None,clusters_0p5=old['greedy_clusters_0p5'],
            archived_std=old['energy_std']))
    with (OUT/'per_molecule.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(per[0]),lineterminator='\n');w.writeheader();w.writerows(per)
    ref={r['mol_id']:r for r in per if r['method']=='chembl3d_gt_pb'}
    summary=[]
    for label in labels:
        rs=[r for r in per if r['method']==label['method'] and r['n_energy']>0]
        std=[r['std_energy'] for r in rs]
        summary.append(dict(method=label['method'],label=label['display_label'],n_molecules=len(rs),n_conformers=sum(r['n_energy'] for r in rs),
            n_unconverged=sum(r['n_unconverged'] for r in rs),mean_clusters_0p5=float(np.mean([r['clusters_0p5'] for r in rs])),
            median_energy_std=float(np.median(std)),mean_energy_std=float(np.mean(std)),median_mean_energy=float(np.median([r['mean_energy'] for r in rs])),
            median_mean_energy_delta_chembl=float(np.median([r['mean_energy']-ref[r['mol_id']]['mean_energy'] for r in rs])),
            median_archived_std=float(np.median([r['archived_std'] for r in rs]))))
    with (OUT/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]),lineterminator='\n');w.writeheader();w.writerows(summary)
    (OUT/'provenance.json').write_text(json.dumps(dict(prepared_on='2026-10-06',protocol=CONFIG,fingerprint=fingerprint,input_files=input_files,
        database=str(DB),mapping=str(DEFAULT_CHEMBL_MAP_CSV),reference_reconstruction='Identity-matched loader; 94/94 core entries have zero archived PB rejections; counts and archived energy summaries verified.'),indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
