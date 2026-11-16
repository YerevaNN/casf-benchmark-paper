"""Render Figure 5 entirely from archived curves; no geometry calculation."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT=Path(__file__).resolve().parent
P=OUT.parents[1]
STYLE=json.loads((P/'sampling-data/palette.json').read_text())
methods=json.loads((P/'sampling-data/display.json').read_text())['methods']
pl_methods=['qwen_1p7b_v2_coordtoken_step41054','loqi_druglike','flowr_druglike','nextmol_dmt_l_druglike']
names=['Qwen','LoQI','FlowR','NExT-Mol']
colors=[STYLE[m.removesuffix('_fixed')][1] for m in methods]
core=pd.read_csv(P/'sampling-data/recovery_curves.csv')
ref=pd.read_csv(P/'results-data/core-ref-protocol/ref/recovery_curves.csv')
budgets=sorted(core.budget.unique()); assert budgets[0]==10 and budgets[-1]==1000
per=pd.read_csv(P/'plinder22-data/per-molecule.csv')
summary=pd.read_csv(OUT/'plinder-curves.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C','xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
fig=plt.figure(figsize=(11.8,8.7)); gs=fig.add_gridspec(2,2,left=.08,right=.85,top=.86,bottom=.19,hspace=.66,wspace=.86)
axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1])]
for ax in axes:
 ax.set_xlim(10,1000); ax.set_xticks([10,100,250,500,1000],['10','100','250','500','1,000'])
 ax.set_ylim(40,100);ax.set_yticks([40,50,60,70,80,90,100]);ax.grid(color='#E7EBEE',lw=.6);ax.set_axisbelow(True)
 ax.set_xlabel('Conformer budget')

def endpoint_labels(ax,values,precision=1,unit='%',gap=3.8):
 ordered=sorted(enumerate(values),key=lambda v:v[1]);pos=[]
 for _,val in ordered:pos.append(max(val,pos[-1]+gap) if pos else val)
 shift=max(0,pos[-1]-(ax.get_ylim()[1]-gap/2))
 for (i,val),loc in zip(ordered,pos):
  ax.plot(1000,val,'o',ms=3.5,color=colors[i],clip_on=False)
  ax.annotate(f'{names[i]}  {val:.{precision}f}{unit}',xy=(1000,val),xytext=(1.025,loc-shift),textcoords=ax.get_yaxis_transform(),ha='left',va='center',fontsize=8,color=colors[i],annotation_clip=False,arrowprops=dict(arrowstyle='-',color=colors[i],lw=.7))
ax=axes[0];ends=[]
for m,col in zip(pl_methods,colors):
 c=summary[summary.method.eq(m)];ax.plot(c.budget,100*c.coverage,lw=1.9,color=col);ends.append(100*c.loc[c.budget.eq(1000),'coverage'].iloc[0])
ax.set_title('A   PLINDER · reference coverage',loc='left',pad=12)
ax.set_ylabel('Expected COV-R (%)');ax.set_xlabel('Supplied conformers sampled');endpoint_labels(ax,ends)
ax=axes[1];ends=[]
for m,col in zip(pl_methods,colors):
 c=summary[summary.method.eq(m)];ax.plot(c.budget,c.mat_r,lw=1.9,color=col);ends.append(c.loc[c.budget.eq(1000),'mat_r'].iloc[0])
ax.set_ylim(.3,1.0);ax.set_yticks([.3,.4,.5,.6,.7,.8,.9,1.0])
ax.set_title('B   PLINDER · matching distance',loc='left',pad=12)
ax.set_ylabel('Expected MAT-R (Å; lower is better)');ax.set_xlabel('Supplied conformers sampled')
endpoint_labels(ax,ends,precision=3,unit=' Å',gap=.05)
counts=pd.read_csv(P/'results-data/core-ref-protocol/chembl-candidate-counts.csv')
for i,(cohort,curves,n) in enumerate([('core',core,94),('ref',ref,1044)],2):
 ax=axes[i];ends=[]
 for m,col in zip(methods,colors):
  c=curves[curves.method.eq(m)&curves.threshold.eq(.75)];ax.plot(c.budget,100*c.recovery,lw=1.9,color=col);ends.append(100*c.loc[c.budget.eq(1000),'recovery'].iloc[0])
 ax.set_title(f'{chr(65+i)}   CASF {cohort} · {n:,} entries',loc='left',pad=30)
 ax.set_ylabel('Expected recovery (%)');ax.set_xlabel('Candidate budget')
 mean=counts[counts.cohort.eq(cohort)].chembl_count.mean()
 ax.axvline(mean,color='#626C78',ls=(0,(3,3)),lw=1.1)
 ax.annotate(f'Mean ChEMBL count: {mean:.1f}',xy=(mean,1.015),xycoords=ax.get_xaxis_transform(),ha='left',va='bottom',fontsize=8,color='#626C78')
 src=P/('sampling-data/reference.json' if cohort=='core' else 'results-data/core-ref-protocol/ref/reference.json')
 reference=json.loads(src.read_text())['hit_rates']['0.75']
 ax.axhline(100*reference,color=STYLE['chembl3d_gt_pb'][1],ls=(0,(4,3)),lw=1.2)
 endpoint_labels(ax,ends)
fig.suptitle('Sampling budget and recovery of observed bound conformations',x=.08,ha='left',y=.975,fontsize=13)
fig.text(.08,.935,'PLINDER: 22 molecules · coverage cutoff: 0.75 Å · MAT-R uses continuous matching distances',fontsize=9,color='#626C78')
handles=[Line2D([],[],color=col,lw=1.9,label=STYLE[m.removesuffix('_fixed')][0]) for m,col in zip(methods,colors)]
handles.append(Line2D([],[],color=STYLE['chembl3d_gt_pb'][1],ls=(0,(4,3)),lw=1.2,label='Stored ChEMBL3D-PB (CASF)'))
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.08),ncol=3,frameon=False,fontsize=8)
fig.text(.08,.048,'PLINDER: supplied ensembles, without a common PoseBusters filter. CASF: rejected candidates consume budget.',fontsize=8,color='#626C78')
fig.text(.08,.025,'Uniform subsets without replacement; sample counts are capped at each available pool. No new conformers generated.',fontsize=8,color='#626C78')
for ext in ['png','pdf','svg']:fig.savefig(P/'figures'/f'figure-budget-combined.{ext}',dpi=240)
plt.close(fig)
svg=P/'figures/figure-budget-combined.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
inputs=[OUT/'plinder-curves.csv', OUT/'plinder-per-molecule-curves.csv',
 P/'sampling-data/recovery_curves.csv', P/'results-data/core-ref-protocol/ref/recovery_curves.csv',
 P/'results-data/core-ref-protocol/chembl-candidate-counts.csv', P/'sampling-data/reference.json',
 P/'results-data/core-ref-protocol/ref/reference.json', P/'sampling-data/palette.json',
 P/'sampling-data/display.json', Path(__file__)]
outputs=[P/'figures'/f'figure-budget-combined.{ext}' for ext in ['pdf','png','svg']]
(OUT/'rendered-files.json').write_text(json.dumps({
 'inputs':{str(path.relative_to(P)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
 'outputs':{str(path.relative_to(P)):hashlib.sha256(path.read_bytes()).hexdigest() for path in outputs}
},indent=2)+'\n')
print('Rendered Figure 5: PLINDER COV-R/MAT-R and CASF core/ref recovery.')
