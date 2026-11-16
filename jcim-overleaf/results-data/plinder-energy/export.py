"""Verify and aggregate saved PLINDER common-hydrogen energies; no rescoring."""
from pathlib import Path
import subprocess,hashlib,json,ast
import pandas as pd,numpy as np
O=Path(__file__).resolve().parent;P=O.parents[1];B=Path('/mnt/weka/mbedrosian/code/casf-benchmark')
source=P/'plinder22-data/per-molecule.csv';d=pd.read_csv(source)
assert len(d)==220 and d.label.nunique()==10
assert d.energy_count.eq(d.n_pb_pass).all()
assert d.energy_std.notna().eq(d.energy_count.gt(0)).all()
cols=['label','smiles','name','n_generated','n_pb_pass','energy_count','energy_min','energy_max','energy_median','energy_std','greedy_clusters_0p5']
d[cols].to_csv(O/'per-molecule.csv',index=False)
style=json.loads((P/'sampling-data/palette.json').read_text())
points=pd.read_csv(P/'plinder22-data/scatter-points.csv').set_index('method')
summary=[]
for method,g in d.groupby('label',sort=False):
 assert len(g)==g.smiles.nunique()==22
 valid=g[g.energy_count.gt(0)]
 summary.append(dict(method=method,label=points.loc[method,'name'],n_molecules=len(valid),n_total=len(g),n_empty=int(g.energy_count.eq(0).sum()),n_conformers=int(g.energy_count.sum()),n_pb_pass=int(g.n_pb_pass.sum()),n_generated=int(g.n_generated.sum()),n_failed=int((g.n_pb_pass-g.energy_count).sum()),mean_clusters_0p5=valid.greedy_clusters_0p5.mean(),median_energy_std=valid.energy_std.median(),mean_energy_std=valid.energy_std.mean()))
pd.DataFrame(summary).to_csv(O/'summary.csv',index=False)
# Match the original evaluation-code hash recorded with the imported PLINDER data.
old=subprocess.check_output(['git','show','85b8fcc74:src/casf_benchmark/analysis/plinder_panel.py'],cwd=B)
original=json.loads((P/'plinder22-data/provenance.json').read_text())['source_sha256']
assert hashlib.sha256(old).hexdigest()==original[str(B/'src/casf_benchmark/analysis/plinder_panel.py')]
(O/'source-plinder-panel.py').write_bytes(old)
functions=['common_hydrogen_mmff_energy','energy_summary','score_conformer_pool']
def get_tree(code,name):return next(n for n in ast.parse(code).body if isinstance(n,ast.FunctionDef) and n.name==name)
current=(B/'src/casf_benchmark/analysis/plinder_panel.py').read_text()
other=Path('/mnt/weka/vtarasov/code/casf-benchmark/src/casf_benchmark/analysis/plinder_panel.py').read_text()
checks={}
for fn in functions:
 expected=ast.dump(get_tree(old.decode(),fn));checks[fn]={'main_matches_original':ast.dump(get_tree(current,fn))==expected,'baseline_checkout_matches_original':ast.dump(get_tree(other,fn))==expected}
assert all(v['main_matches_original'] for v in checks.values())
# Confirm current upstream energy aggregates for all non-RDKit methods.
upstreams={};errors=[]
for method,g in d.groupby('label'):
 if method.startswith('rdkit_'):continue
 f=(Path('/mnt/weka/mbedrosian/pharma_generation_analysis/v2_1p7b_coordtoken_step41054/plinder_panel') if method.startswith('qwen_') else B/'data/results/plinder_panel')/method/'per_molecule.csv'
 src=pd.read_csv(f).set_index('smiles').loc[g.smiles]
 np.testing.assert_allclose(g[['n_pb_pass','energy_count','energy_std']].to_numpy(),src[['n_pb_pass','energy_count','energy_std']].to_numpy(),rtol=1e-12,atol=1e-12)
 upstreams[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
# RDKit raw uses the archived atom-order-corrected validity/energy records.
f=P/'plinder22-data/minimized-runs/rdkit_random_raw_druglike-corrected-validity-energy.csv';src=pd.read_csv(f).set_index('smiles');g=d[d.label.eq('rdkit_random_raw_druglike')]
np.testing.assert_allclose(g[['n_pb_pass','energy_count','energy_std']],src.loc[g.smiles,['n_pb_pass','energy_count','energy_std']],rtol=1e-12,atol=1e-12)
upstreams[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
f=Path('/mnt/weka/vtarasov/code/casf-benchmark/data/results/plinder_panel/rdkit_random_minimized_22_druglike/per_molecule.csv');src=pd.read_csv(f).set_index('smiles');g=d[d.label.eq('rdkit_random_minimized_22_druglike')]
np.testing.assert_allclose(g[['n_pb_pass','energy_count','energy_std']],src.loc[g.smiles,['n_pb_pass','energy_count','energy_std']],rtol=1e-12,atol=1e-12)
upstreams[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
protocol=json.loads((P/'energy-data/protocol.json').read_text());protocol={k:v for k,v in protocol.items() if k not in ['workers','ligand_set','target','source_database_sha256']}
info={'protocol':protocol,'energy_pool':'Saved PoseBusters-passing subset; no new PB replay or energy evaluation.','geometry_pool':'Full supplied ensembles for molecules with defined energy SD (21 for FlowR, 22 for others); not the same conformer subset as energy.','aggregation':'Median of per-molecule population energy SD; mean cluster count over the same molecules.','finite_energies':int(d.energy_count.sum()),'pb_passing':int(d.n_pb_pass.sum()),'supplied_conformers':int(d.n_generated.sum()),'source_sha256':{str(source):hashlib.sha256(source.read_bytes()).hexdigest(),**upstreams},'original_evaluator_sha256':hashlib.sha256(old).hexdigest(),'function_comparisons':checks,'figure':'Main Figure 2; PLINDER added above unchanged CASF core/ref measurements.'}
(O/'protocol-audit.json').write_text(json.dumps(info,indent=2)+'\n')
print(pd.DataFrame(summary)[['label','median_energy_std','mean_energy_std']].to_string(index=False));print(checks);print('PLINDER finite / PB-passing:',info['finite_energies'],info['pb_passing'])
