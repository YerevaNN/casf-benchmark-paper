"""Update S1, S2 and S3 from the current archived ref data."""
from pathlib import Path
import pandas as pd
P=Path(__file__).resolve().parents[2]
table3=pd.read_csv(P/'results-data/three-dataset-tables/table3.csv')
r=table3[table3.cohort.eq('ref')].set_index(['method','budget'])
p=P/'tables/supporting-recovery.tex';lines=p.read_text().splitlines()
for i,line in enumerate(lines):
 name=line.split(' & ')[0]
 if name not in r.index.get_level_values('method'):continue
 a=r.loc[(name,'chembl_count')]
 vals=[f'{a.recovery:.1f}',f'{a.mean_retained:.1f}']
 if name=='ChEMBL3D-PB':vals+=['---','---']
 else:
  b=r.loc[(name,'fixed')];vals +=[f'{b.recovery:.1f}',f'{b.mean_retained:.1f}']
 lines[i]=' & '.join([name]+vals)+r' \\'
lines[0]='% Rebuild with results-data/ref-panel/build_supporting.py.'
p.write_text('\n'.join(lines)+'\n')
size=pd.read_csv(P/'results-data/size-flexibility/casf-size-source.csv').set_index(['set','method','heavy_atoms'])
methods=['chembl3d_gt_pb','loqi_raw_fixed','flowr_raw_fixed','qwen_1p7b_v2_coordtoken_step41054_fixed']
groups=[('<20',r'$<20$'),('20-29','20--29'),('30-39','30--39'),('40+',r'$\geq40$')]
p=P/'tables/size.tex';lines=p.read_text().splitlines()
for cohort in ['core','ref']:
 rows=[]
 for group,label in groups:
  d=[size.loc[(cohort,m,group)] for m in methods]
  assert len({v.n for v in d})==1
  rows.append(' & '.join([cohort,label,str(int(d[0].n))]+[f'{100*v.hit_075:.1f}' for v in d])+r' \\')
 position=next(i for i,l in enumerate(lines) if l.startswith(cohort+' &'))
 lines=[l for l in lines if not l.startswith(cohort+' &')];lines[position:position]=rows
lines[0]='% Rebuild with results-data/ref-panel/build_supporting.py.'
p.write_text('\n'.join(lines)+'\n')
# Current primary Qwen values in Supporting Table S1.
tiers=pd.read_csv(P/'figure-data/recovery-table-tiers.csv')
d=tiers[tiers.cohort.eq('core')&tiers.key.eq('qwen_1p7b_v2_coordtoken_step41054')].set_index('tier')
p=P/'tables/supporting-rmsd.tex';lines=p.read_text().splitlines()
for i,line in enumerate(lines):
 if line.startswith('Qwen 1.7B FSQ &'):
  values=[d.loc[t,c] for t in ['chembl_count','fixed'] for c in ['mean_minimum_rmsd','mean_ensemble_median_rmsd']]
  lines[i]='Qwen 1.7B FSQ & '+' & '.join(f'{v:.3f}' for v in values)+r' \\'
p.write_text('\n'.join(lines)+'\n')
