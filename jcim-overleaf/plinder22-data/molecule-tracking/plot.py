"""Connect the same PLINDER molecule across methods using existing matching scores."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.colors import qualitative

OUT = Path(__file__).resolve().parent
SOURCE = OUT/'molecule_matching_values.csv'
df = pd.read_csv(SOURCE)
methods = json.loads((OUT/'methods.json').read_text())
assert len(df)==198 and not df.duplicated(['label','smiles']).any()
assert df.groupby('name').smiles.nunique().eq(1).all()
assert np.isfinite(df[['mat_r','mat_p']]).all().all()
assert df.groupby('smiles').label.nunique().eq(9).all()
wide = {metric:df.pivot(index='name',columns='label',values=metric).loc[:,list(methods)] for metric in ['mat_r','mat_p']}
outliers = pd.concat([df.loc[df.groupby('label')[metric].idxmax(),['label','name',metric]].rename(columns={metric:'value'}).assign(metric=metric) for metric in wide],ignore_index=True)
outliers.to_csv(OUT/'largest_outliers.csv',index=False)
focus = sorted(set(outliers.name))
assert len(focus)==3
focus_colors = {'Actinonin':'#0072B2','Amprenavir':'#CE6428','Imatinib':'#AA589B'}
all_names = focus + sorted(set(df.name)-set(focus))
other_colors = [c for c in qualitative.Dark24 if c not in focus_colors.values()]
colors = {**focus_colors, **dict(zip(all_names[len(focus):], other_colors))}

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C','xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
fig,axes=plt.subplots(2,1,figsize=(13.5,9.5),sharex=True)
fig.subplots_adjust(left=.09,right=.84,top=.85,bottom=.18,hspace=.28)
x=np.arange(len(methods))
last=len(methods)-1
for ax,metric,title,upper in zip(axes,wide,['A   Experimental → generated: MAT-R','B   Generated → experimental: MAT-P'],[2.55,3.85]):
    values=wide[metric]
    for name in all_names:
        if name in focus:continue
        ax.plot(x,values.loc[name],color='#CED3D8',lw=.8,marker='o',ms=2.4,zorder=1)
    for name in focus:
        vals=values.loc[name].to_numpy()
        ax.plot(x,vals,color=colors[name],lw=1.8,marker='o',ms=4,zorder=3)
    # Label the worst molecule for each method with its measured value.
    for j,method in enumerate(methods):
        name=values[method].idxmax();v=values.loc[name,method]
        ax.scatter([j],[v],s=86,facecolors='none',edgecolors=colors[name],lw=1.2,zorder=4)
        ax.annotate(f'{v:.2f}',(j,v),xytext=(0,10),textcoords='offset points',ha='center',fontsize=8,color=colors[name])
    # End labels for highlighted trajectories, spaced to avoid collisions.
    ends=sorted([(values.loc[name].iloc[-1],name) for name in focus])
    spaced=[]
    for val,name in ends:
        pos=max(val,spaced[-1][0]+upper*.065) if spaced else val
        spaced.append((pos,name,val))
    for pos,name,val in spaced:
        ax.annotate(name,(last,val),xytext=(last+.28,pos),textcoords='data',ha='left',va='center',color=colors[name],fontsize=9,annotation_clip=False,arrowprops={'arrowstyle':'-','color':colors[name],'lw':.65})
    ax.set(xlim=(-.15,last+.15),ylim=(0,upper),ylabel='RMSD (Å; lower is better)')
    ax.set_title(title,loc='left',fontsize=11,pad=10)
    ax.grid(axis='y',color='#E7EBEE',lw=.6);ax.set_axisbelow(True)
axes[1].set_xticks(x,list(methods.values()))
for tick,method in zip(axes[1].get_xticklabels(),methods):
    if 'casf_dedup' in method: tick.set_fontweight('bold')
fig.suptitle('Tracking individual molecules across conformer generators',x=.09,y=.975,ha='left',fontsize=14)
fig.text(.09,.935,'PLINDER-22 · one connected line per molecule · the same molecule keeps its color in both panels',fontsize=10,color='#626C78')
fig.text(.09,.90,'Three molecules are worst for at least one method (colored); the other 19 are grey. New Qwen run is bold.',fontsize=9,color='#626C78')
fig.text(.09,.087,'Outlined markers identify each method’s largest distance; numbers give that molecule’s value in Å.',fontsize=9,color='#626C78')
fig.text(.09,.052,'Full supplied ensembles; no common PoseBusters filter. † Torsion raw starts from an experimental conformation.',fontsize=8,color='#626C78')
fig.text(.09,.022,'Methods are categorical, ordered by mean MAT-R. Connecting lines track identity, not a continuous trend.',fontsize=8,color='#626C78')
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'plinder22-molecule-tracking.{ext}',dpi=250)
plt.close(fig)

# Standalone interactive companion: any molecule can be isolated in both panels.
labels=[v.replace('\n',' ') for v in methods.values()]
p=make_subplots(rows=2,cols=1,shared_xaxes=True,vertical_spacing=.12,subplot_titles=['Experimental → generated: MAT-R','Generated → experimental: MAT-P'])
for i,metric in enumerate(wide,1):
    for name in all_names:
        p.add_trace(go.Scatter(x=labels,y=wide[metric].loc[name],name=name,legendgroup=name,showlegend=i==1,mode='lines+markers',line={'color':colors[name],'width':2 if name in focus else 1},marker={'size':6},hovertemplate=f'<b>{name}</b><br>%{{x}}<br>{metric.upper().replace("_","-")}: %{{y:.3f}} Å<extra></extra>'),row=i,col=1)
    p.update_yaxes(title_text=metric.upper().replace('_','-')+' (Å; lower is better)',rangemode='tozero',row=i,col=1)
visible_focus=[True if name in focus else 'legendonly' for _ in wide for name in all_names]
p.update_layout(title={'text':'Track each PLINDER-22 molecule across methods<br><sup>Hover for values. Click a legend name to toggle it; double-click to isolate it in both panels.</sup>','x':.06},template='plotly_white',height=980,width=1450,margin={'t':165,'b':120,'r':240,'l':85},legend={'title_text':'Molecule','groupclick':'togglegroup','tracegroupgap':0},updatemenus=[{'type':'buttons','direction':'right','x':0,'y':1.14,'buttons':[{'label':'Show all 22','method':'restyle','args':[{'visible':True}]},{'label':'Highlight outlier molecules','method':'restyle','args':[{'visible':visible_focus}]}]}])
p.add_annotation(text='Full supplied ensembles; no common PoseBusters filter. † Torsion raw starts from an experimental conformation.',xref='paper',yref='paper',x=0,y=-.16,showarrow=False,font={'size':11},xanchor='left')
p.write_html(OUT/'plinder22-molecule-tracking.html',include_plotlyjs=True,full_html=True)
(OUT/'caption.md').write_text('Molecule-level matching distances on the cleaned PLINDER 22-molecule panel (1,067 reference conformers). Each line connects values for the same molecule across nine methods: the eight previously displayed methods plus Qwen 0.6B CASF-dedup binned PM12 step64000. The original Qwen 1.7B FSQ step47023 remains separately identified. Panels report MAT-R and MAT-P, respectively; lower is better. Methods are ordered by overall mean MAT-R and are categorical. In the static figure, the union of each method’s worst-performing molecule in either metric is highlighted in color (Actinonin, Amprenavir, and Imatinib), while the other 19 are grey. Outlined points mark the largest distance for each method and metric, with numerical labels in Å. The interactive companion assigns every molecule a distinct color and supports linked legend selection across panels. Existing supplied-ensemble measurements are used without new generation, optimization, or RMSD calculations. There is no common PoseBusters filter. The torsion baseline starts from an experimental reference. The CASF-dedup label identifies the new run; this plot does not audit its training exclusions. Bortezomib is excluded from the cleaned panel.\n')
(OUT/'provenance.json').write_text(json.dumps({'source':str(SOURCE),'sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'highlight_rule':'union of per-method maxima over MAT-R and MAT-P','highlighted_molecules':focus,'structural_recalculation':False,'paper_updated':False},indent=2)+'\n')
print(outliers.to_string(index=False))
