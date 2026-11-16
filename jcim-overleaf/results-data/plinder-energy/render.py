"""Render Figure 2 with PLINDER above the unchanged CASF energy panels."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter,NullLocator

HERE=Path(__file__).resolve().parent
PAPER=HERE.parents[1]
style=json.loads((PAPER/'sampling-data/palette.json').read_text())
colors={name:color for name,color,_ in style.values()}
core=pd.read_csv(PAPER/'energy-data/summary.csv')
ref=pd.read_csv(PAPER/'results-data/core-ref-protocol/ref/summary.csv')
plinder=pd.read_csv(HERE/'summary.csv')
order=core.label.tolist()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':10.5,
 'axes.labelsize':9.5,'xtick.labelsize':8.5,'ytick.labelsize':9,
 'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C',
 'xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none',
 'savefig.facecolor':'white','figure.facecolor':'white'})
fig,axes=plt.subplots(3,2,figsize=(9.5,10.2))
fig.subplots_adjust(left=.235,right=.97,top=.90,bottom=.11,hspace=.58,wspace=.33)
points=[]
for i,(cohort,data,title) in enumerate([
 ('plinder',plinder,'PLINDER · supplied ensembles'),
 ('core',core,'CASF core · C-count target'),('ref',ref,'CASF ref · C-count target')]):
 methods=[m for m in order if m in set(data.label)]
 for j,metric in enumerate(['mean_clusters_0p5','median_energy_std']):
  ax=axes[i,j];ax.set_ylim(len(methods)-.4,-.6)
  ax.set_yticks(range(len(methods)),methods if j==0 else []);ax.tick_params(axis='y',length=0)
  ax.grid(axis='x',color='#E7EBEE',lw=.6);ax.set_axisbelow(True)
  if cohort!='plinder':ax.axhline(len(methods)-1.5,color='#7B838D',lw=.7,ls=(0,(4,3)))
  panel_title='PLINDER · PB-passing subsets' if cohort=='plinder' and j==1 else title
  ax.set_title(f'{chr(65+2*i+j)}   {panel_title}',loc='left',pad=12)
  if j==0:
   ax.set_xlim(0,500 if cohort=='plinder' else 66)
   ax.set_xticks([0,100,200,300,400,500] if cohort=='plinder' else [0,20,40,60])
   ax.set_xlabel('Mean clusters at 0.5 Å',labelpad=6)
  elif cohort=='plinder':
   ax.set_xscale('log');ax.set_xlim(1,10000);ax.set_xticks([1,10,100,1000,10000])
   ax.xaxis.set_major_formatter(ScalarFormatter());ax.xaxis.set_minor_locator(NullLocator())
   ax.set_xlabel('Median energy SD (kcal/mol; log scale)',labelpad=6)
  else:
   ax.set_xlim(0,49);ax.set_xticks([0,10,20,30,40])
   ax.set_xlabel('Median energy SD (kcal/mol)',labelpad=6)
  for k,m in enumerate(methods):
   r=data[data.label.eq(m)].iloc[0];value=float(r[metric]);left=1 if cohort=='plinder' and j==1 else 0
   ax.barh(k,value-left,left=left,color=colors[m],height=.58,zorder=3)
   offset=value*1.13 if left else value+(.9 if cohort!='plinder' else 7)
   ax.text(offset,k,f'{value:,.1f}',va='center',fontsize=8.5)
   points.append(dict(cohort=cohort,method=m,metric=metric,value=value,n_molecules=int(r.n_molecules)))
fig.suptitle('Geometric diversity and energy dispersion',x=.235,ha='left',y=.978,fontsize=13)
fig.text(.235,.944,'Common hydrogen preparation · MMFF94s · no coordinate optimization',fontsize=9,color='#626C78')
fig.text(.235,.055,'PLINDER clusters: full supplied pools; energies: PoseBusters-passing subsets.',fontsize=8.5,color='#626C78')
fig.text(.235,.032,'Paired molecule sets: 21 for FlowR, 22 for other PLINDER methods. Panel axes differ.',fontsize=8.5,color='#626C78')
outputs=[]
for ext in ['pdf','png','svg']:
 path=PAPER/'figures'/f'figure-energy-dispersion.{ext}';fig.savefig(path,dpi=260);outputs.append(path)
plt.close(fig)
svg=outputs[-1];svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
pd.DataFrame(points).to_csv(HERE/'figure-points.csv',index=False)
sources=[HERE/'summary.csv',HERE/'per-molecule.csv',HERE/'protocol-audit.json',
 PAPER/'energy-data/summary.csv',PAPER/'results-data/core-ref-protocol/ref/summary.csv',
 PAPER/'sampling-data/palette.json',Path(__file__)]
(HERE/'rendered-files.json').write_text(json.dumps({'inputs':{str(p.relative_to(PAPER)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'outputs':{str(p.relative_to(PAPER)):hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs}},indent=2)+'\n')
print('Rendered Figure 2: PLINDER and CASF energy/diversity.')
