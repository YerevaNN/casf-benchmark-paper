"""Render the selected sampling methods in one two-panel figure."""
from pathlib import Path
import json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT=Path(__file__).resolve().parent
STYLE=json.loads((OUT/'palette.json').read_text())
DISPLAY=json.loads((OUT/'display.json').read_text())
curves=pd.read_csv(OUT/'recovery_curves.csv')
frontier=pd.read_csv(OUT/'recovery80_frontier.csv')
ref=json.loads((OUT/'reference.json').read_text())
methods=DISPLAY['methods'];threshold=DISPLAY['hit_threshold']
assert set(methods)<=set(curves.method) and set(methods)<=set(frontier.method)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.titleweight':'medium','axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C',
 'xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
 'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,4.6))
fig.subplots_adjust(left=.09,right=.97,top=.80,bottom=.29,wspace=.43)
for ax in (a,b):
    ax.set_xscale('log');ax.set_xlim(10,1000)
    ax.set_xticks([10,100,1000],labels=['10','100','1,000'])
    ax.minorticks_off();ax.grid(color='#E7EBEE',lw=.6);ax.set_axisbelow(True)
    ax.set_xlabel('Candidate budget')
for method in methods:
    key=method.removesuffix('_fixed');label,color,_=STYLE[key]
    c=curves[curves.method.eq(method)&curves.threshold.eq(threshold)]
    f=frontier[frontier.method.eq(method)]
    assert len(c)==len(f)>0
    a.fill_between(c.budget,100*c.ci_low,100*c.ci_high,color=color,alpha=.09,lw=0)
    a.plot(c.budget,100*c.recovery,color=color,lw=1.9)
    b.plot(f.budget,f.rmsd_threshold,color=color,lw=1.9)
reference_color=STYLE['chembl3d_gt_pb'][1]
a.axhline(100*ref['hit_rates'][str(threshold)],color=reference_color,ls=(0,(4,3)),lw=1.2)
b.axhline(ref['rmsd_threshold'],color=reference_color,ls=(0,(4,3)),lw=1.2)
a.set_title(f'A   Recovery at {threshold:g} Å',loc='left',pad=12)
b.set_title('B   Threshold for 80% recovery',loc='left',pad=12)
a.set_ylabel('Expected recovery (%)');a.set_ylim(0,100);a.set_yticks([0,20,40,60,80,100])
b.set_ylabel('RMSD threshold (Å)');b.set_ylim(0,1.5);b.set_yticks([0,.3,.6,.9,1.2,1.5])
fig.suptitle('Sampling budget, recovery, and matching tolerance',x=.09,ha='left',y=.96,fontsize=11.5)
fig.text(.09,.89,'94 CASF entries · saved candidate pools · higher in A and lower in B are better',fontsize=8,color='#626C78')
handles=[Line2D([],[],color=STYLE[m.removesuffix('_fixed')][1],lw=1.9,label=STYLE[m.removesuffix('_fixed')][0]) for m in methods]
handles.append(Line2D([],[],color=reference_color,ls=(0,(4,3)),lw=1.2,label='Stored ChEMBL3D-PB'))
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.07),ncol=3,frameon=False,fontsize=8,columnspacing=1.5,handlelength=2.7)
fig.text(.5,.025,'Shading in A: pointwise 95% entry-bootstrap intervals. Rejected candidates consume budget.',ha='center',fontsize=7.2,color='#626C78')
for ext in ['pdf','svg','png']:fig.savefig(OUT/(DISPLAY['figure_name']+'.'+ext),dpi=350)
svg=OUT/(DISPLAY['figure_name']+'.svg');svg.write_text('\n'.join(s.rstrip() for s in svg.read_text().splitlines())+'\n')
plt.close(fig)
print('Rendered the combined four-method sampling figure.')
