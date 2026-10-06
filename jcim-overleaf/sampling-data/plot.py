"""Render the two review alternatives from the archived finite-pool analysis."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT=Path(__file__).resolve().parent
STYLE=json.loads((OUT/'palette.json').read_text())
curves=pd.read_csv(OUT/'recovery_curves.csv');frontier=pd.read_csv(OUT/'recovery80_frontier.csv');ref=json.loads((OUT/'reference.json').read_text())
methods=list(curves.method.unique())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.titleweight':'medium','axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C','xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})
def style(method):
    key=method.removesuffix('_fixed');label,color,marker=STYLE[key]
    return label,color,'--' if key.endswith('_minimized') else '-'
def axis(ax):
    ax.set_xscale('log');ax.set_xlim(10,1000);ax.set_xticks([10,100,1000],labels=['10','100','1,000'])
    ax.minorticks_off();ax.grid(color='#E7EBEE',lw=.6);ax.set_axisbelow(True)
    ax.set_xlabel('Candidate budget')
def legend(fig,y):
    handles=[Line2D([],[],color=style(m)[1],ls=style(m)[2],lw=1.7,label=style(m)[0]) for m in methods]
    handles.append(Line2D([],[],color=STYLE['chembl3d_gt_pb'][1],ls=(0,(4,3)),lw=1.3,label='Stored ChEMBL3D-PB'))
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,y),ncol=3,frameon=False,fontsize=7.5,columnspacing=1.5,handlelength=2.7)
def save(fig,name):
    for ext in ['pdf','svg','png']:fig.savefig(OUT/(name+'.'+ext),dpi=350)
    svg=OUT/(name+'.svg');svg.write_text('\n'.join(s.rstrip() for s in svg.read_text().splitlines())+'\n')
    plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(7.2,5.2),sharey=True)
fig.subplots_adjust(left=.09,right=.98,top=.80,bottom=.34,wspace=.27)
for ax,t,letter in zip(axes,[.5,.75,1.],['A','B','C']):
    axis(ax);ax.set_title(f'{letter}   Hit@{t:g} Å',loc='left',pad=12)
    for method in methods:
        d=curves[curves.method.eq(method)&curves.threshold.eq(t)]
        label,color,ls=style(method)
        ax.fill_between(d.budget,100*d.ci_low,100*d.ci_high,color=color,alpha=.045,lw=0)
        ax.plot(d.budget,100*d.recovery,color=color,ls=ls,lw=1.65)
    ax.axhline(100*ref['hit_rates'][str(t)],color=STYLE['chembl3d_gt_pb'][1],ls=(0,(4,3)),lw=1.2)
    ax.set_ylim(0,100);ax.set_yticks([0,20,40,60,80,100])
axes[0].set_ylabel('Expected recovery (%)')
fig.suptitle('Sampling budget and matching tolerance',x=.09,ha='left',y=.96,fontsize=12)
fig.text(.09,.89,'Uniform subsets of the archived candidate pools · 94 CASF entries',fontsize=8,color='#626C78')
legend(fig,.085)
fig.text(.5,.03,'Shading: pointwise 95% molecule-bootstrap intervals; rejected candidates consume budget.',ha='center',fontsize=7.2,color='#626C78')
save(fig,'figure-budget-thresholds')
fig,ax=plt.subplots(figsize=(7.2,5.3));fig.subplots_adjust(left=.12,right=.97,top=.84,bottom=.34)
axis(ax)
for method in methods:
    d=frontier[frontier.method.eq(method)];label,color,ls=style(method)
    ax.plot(d.budget,d.rmsd_threshold,color=color,ls=ls,lw=1.75)
ax.axhline(ref['rmsd_threshold'],color=STYLE['chembl3d_gt_pb'][1],ls=(0,(4,3)),lw=1.3)
ax.set_ylabel('RMSD threshold for 80% recovery (Å)')
ax.set_ylim(0,max(frontier.rmsd_threshold.max(),ref['rmsd_threshold'])*1.08)
fig.suptitle('How closely can each budget recover 80% of molecules?',x=.12,ha='left',y=.96,fontsize=11)
fig.text(.12,.89,'Lower curves indicate closer matches at the same candidate budget',fontsize=8,color='#626C78')
legend(fig,.075)
fig.text(.5,.025,'Exact threshold for 80% expected recovery within the saved pools; no extrapolation beyond 1,000.',ha='center',fontsize=7.1,color='#626C78')
save(fig,'figure-budget-recovery80')
print('Rendered both sampling alternatives.')
