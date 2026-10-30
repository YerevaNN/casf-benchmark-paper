"""Update S2 and the ref rows of S3 from the current archived ref data."""
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
size=pd.read_csv(P/'results-data/size-flexibility/casf-size-source.csv');size=size[size['set'].eq('ref')].set_index(['method','heavy_atoms'])
methods=['chembl3d_gt_pb','loqi_raw_fixed','flowr_raw_fixed','qwen_1p7b_fsq_bigdata_step47023_fixed']
groups=[('<20',r'$<20$'),('20-29','20--29'),('30-39','30--39'),('40+',r'$\geq40$')]
rows=[]
for group,label in groups:
 d=[size.loc[(m,group)] for m in methods]
 assert len({v.n for v in d})==1
 rows.append(' & '.join(['ref',label,str(int(d[0].n))]+[f'{100*v.hit_075:.1f}' for v in d])+r' \\')
p=P/'tables/size.tex';lines=p.read_text().splitlines();position=next(i for i,l in enumerate(lines) if l.startswith('ref &'))
lines=[l for l in lines if not l.startswith('ref &')];lines[position:position]=rows
lines[0]='% Core: original size analysis. Ref: results-data/ref-panel/build_supporting.py.'
p.write_text('\n'.join(lines)+'\n')
