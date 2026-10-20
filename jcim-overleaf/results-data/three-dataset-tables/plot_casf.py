"""CASF core/ref scatter using the archived table data and manuscript style."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[2]
OUT = PAPER / 'figures'
DATA = PAPER / 'results-data/three-dataset-tables'
SOURCES = [DATA/'table1.csv', DATA/'table2.csv', PAPER/'sampling-data/palette.json',
           PAPER/'results-data/diversity-rmsd-points.csv']
palette = json.loads(SOURCES[2].read_text())
styles = {name: (color, marker) for name, color, marker in palette.values()}
diversity = pd.read_csv(SOURCES[0])
recovery = pd.read_csv(SOURCES[1])
points = diversity[['method','cohort','clusters_0p5','n_defined']].merge(
    recovery[['method','cohort','distance','n_total','n_distance']],
    on=['method','cohort'], validate='one_to_one').dropna(subset=['clusters_0p5','distance'])
points = points.loc[points.cohort.isin(['core', 'ref'])].copy()
assert points.groupby('cohort').size().to_dict() == {'core':11,'ref':11}
assert points.n_defined.eq(points.n_distance).all()
assert np.isfinite(points[['clusters_0p5','distance']]).all().all()
old = pd.read_csv(SOURCES[3]).set_index('label')
for r in points.loc[points.cohort.eq('core')].itertuples():
    np.testing.assert_allclose([r.clusters_0p5,r.distance],old.loc[r.method,['clusters_0p5','best_rmsd']].astype(float))
points['distance_metric'] = 'Best RMSD'
points['color'] = points.method.map(lambda m: styles[m][0])
points['marker'] = points.method.map(lambda m: styles[m][1])
points.to_csv(DATA/'casf-scatter-points.csv',index=False)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
    'axes.titleweight':'medium','axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
    'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C',
    'xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,
    'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})

# Point offsets preserve direct colored labels and short leaders from the current figure.
offsets = {
    'core': {
        'ChEMBL3D-PB':(10,-20), 'RDKit raw':(10,3), 'RDKit minimized':(12,4),
        'Torsion raw':(-12,8), 'Torsion minimized':(10,8), 'LoQI':(-10,13),
        'Torsional Diffusion':(10,7), 'MCF drugs-L':(12,9),
        'NExT-Mol DMT-L':(12,-9), 'FlowR':(-12,-4), 'Qwen 1.7B FSQ':(-12,-17)},
    'ref': {
        'ChEMBL3D-PB':(10,-20), 'RDKit raw':(10,-24), 'RDKit minimized':(10,-12),
        'Torsion raw':(10,7), 'Torsion minimized':(10,6), 'LoQI':(10,5),
        'Torsional Diffusion':(16,8), 'MCF drugs-L':(10,12),
        'NExT-Mol DMT-L':(-10,-13), 'FlowR':(12,8), 'Qwen 1.7B FSQ':(10,-17)},
}
fig, axes = plt.subplots(1,2,figsize=(12,5.5))
fig.subplots_adjust(left=.075,right=.975,top=.79,bottom=.22,wspace=.27)
fig.suptitle('Geometric diversity and proximity to observed bound conformations',
             x=.055,y=.975,ha='left',fontsize=13)
fig.text(.055,.928,'CASF core and ref · ChEMBL-count target · PoseBusters-filtered ensembles',
         fontsize=9,color='#626C78')
annotations=[]
for ax, cohort, title, subtitle in zip(axes, ['core','ref'],
        ['A   CASF core · 94 entries','B   CASF ref · 1,236 entries'],
        ['','']):
    ax.set_title(title,loc='left',pad=29,fontsize=11)
    ax.text(0,1.035,subtitle,transform=ax.transAxes,ha='left',fontsize=8,color='#626C78')
    ax.grid(True,color='#E7EBEE',linewidth=.6,zorder=0);ax.set_axisbelow(True)
    ax.tick_params(length=3,width=.6)
    d=points.loc[points.cohort.eq(cohort)]
    for r in d.itertuples():
        color,marker=styles[r.method]
        ax.scatter(r.clusters_0p5,r.distance,s=65,marker=marker,color=color,
                   edgecolor='white',linewidth=.65,zorder=5)
        dx,dy=offsets[cohort][r.method]
        name=r.method
        annotation=ax.annotate(name,(r.clusters_0p5,r.distance),xytext=(dx,dy),
            textcoords='offset points',fontsize=8,ha='right' if dx<0 else 'left',
            va='bottom',color=color,
            arrowprops={'arrowstyle':'-','color':color,'lw':.55,'alpha':.6})
        annotations.append((cohort,r.method,annotation))
    ax.set(xlim=(15,62),ylim=(.45,.74),ylabel='Best RMSD (Å; lower is better)')
    ax.set_xticks([20,30,40,50,60]);ax.set_yticks([.45,.50,.55,.60,.65,.70])
    baseline=d.loc[d.method.eq('ChEMBL3D-PB')].iloc[0]
    ax.axvline(baseline.clusters_0p5,color='#7B838D',lw=.85,ls=(0,(4,3)),zorder=1)
    ax.axhline(baseline.distance,color='#7B838D',lw=.85,ls=(0,(4,3)),zorder=1)
    ax.set_xlabel('Mean geometric clusters at 0.5 Å',labelpad=9)

fig.text(.52,.065,'Dashed guides: stored ChEMBL3D-PB within each CASF panel · lower right: more clusters and lower matching distance',
         ha='center',fontsize=8,color='#626C78')
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'figure-casf-diversity-recovery.{ext}',dpi=300)
plt.close(fig)

# Keep generated vector files clean in source control.
svg = PAPER / 'figures/figure-casf-diversity-recovery.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
