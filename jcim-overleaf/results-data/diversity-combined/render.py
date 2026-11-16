"""Render main Figure 1 and Supporting Figure S8 from current archived measurements."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
PAPER=HERE.parents[1]
OUT=PAPER/'figures'
PL=PAPER/'plinder22-data'
DATA=PAPER/'results-data/three-dataset-tables'
styles=json.loads((PAPER/'sampling-data/palette.json').read_text())
methods = {
    'qwen_1p7b_v2_coordtoken_step41054':'qwen_1p7b_v2_coordtoken_step41054',
    'loqi_druglike':'loqi_raw', 'torsional_diffusion_druglike':'torsional_diffusion_raw',
    'mcf_drugs_l_druglike':'mcf_drugs_l_raw', 'flowr_druglike':'flowr_raw',
    'nextmol_dmt_l_druglike':'nextmol_dmt_l_raw',
    'rdkit_random_raw_druglike':'rdkit_random_raw', 'torsion_independent_raw_druglike':'torsion_raw',
    'rdkit_random_minimized_22_druglike':'rdkit_random_minimized',
    'torsion_independent_minimized_druglike':'torsion_minimized',
}
curve_methods = [
    'rdkit_random_raw_druglike', 'mcf_drugs_l_druglike', 'flowr_druglike',
    'torsional_diffusion_druglike', 'nextmol_dmt_l_druglike', 'loqi_druglike',
    'qwen_1p7b_v2_coordtoken_step41054',
]
summary = pd.read_csv(PL / 'summary.csv').set_index('label')
points = pd.DataFrame([dict(method=m, name=styles[k][0], color=styles[k][1],
    marker=styles[k][2], clusters_0p5=summary.loc[m, 'mean_clusters_0p5'],
    mat_r=summary.loc[m, 'mat_r_mean']) for m,k in methods.items()])
points.to_csv(HERE / 'plinder-points.csv', index=False)
curves = pd.read_csv(PL / 'threshold-curves.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,
    'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
    'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C',
    'xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,
    'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})
plinder_points=points.copy()
diversity=pd.read_csv(DATA/'table1.csv')
recovery=pd.read_csv(DATA/'table2.csv')
points=diversity[['method','cohort','clusters_0p5','n_defined']].merge(
 recovery[['method','cohort','distance','n_total','n_distance']],on=['method','cohort'],validate='one_to_one')
points=points[points.cohort.isin(['core','ref'])].copy()
assert points.groupby('cohort').size().to_dict()=={'core':11,'ref':11}
assert points.n_defined.eq(points.n_distance).all()
assert np.isfinite(points[['clusters_0p5','distance']]).all().all()
points.to_csv(HERE/'casf-points.csv',index=False)
fig=plt.figure(figsize=(12,9.4))
grid=fig.add_gridspec(2,2,left=.08,right=.97,top=.88,bottom=.13,hspace=.63,wspace=.30)
axes=[fig.add_subplot(grid[0,:]),fig.add_subplot(grid[1,0]),fig.add_subplot(grid[1,1])]
fig.suptitle('Geometric diversity and proximity to observed bound conformations',x=.08,y=.975,ha='left',fontsize=13)
for ax in axes:
 ax.grid(color='#E7EBEE',lw=.6);ax.set_axisbelow(True);ax.tick_params(length=3,width=.6)
a=axes[0]
a.set_title('A   PLINDER · 22 molecules · supplied ensembles',loc='left',pad=14)
offsets={'RDKit minimized':(-8,7),'Torsion minimized':(10,-12),'RDKit raw':(10,7),'Torsion raw':(-10,12),'LoQI':(10,8),
         'Torsional Diffusion':(10,-15),'MCF drugs-L':(10,17),
         'NExT-Mol DMT-L':(10,5),'FlowR':(-12,7),'Qwen 1.7B FSQ':(-12,-10)}
for r in plinder_points.itertuples():
    a.scatter(r.clusters_0p5,r.mat_r,s=65,marker=r.marker,color=r.color,
              edgecolor='white',linewidth=.65,zorder=5)
    dx,dy=offsets[r.name]
    a.annotate(r.name,
               (r.clusters_0p5,r.mat_r),xytext=(dx,dy),textcoords='offset points',
               ha='right' if dx<0 else 'left',va='bottom',fontsize=8,color=r.color,
               arrowprops={'arrowstyle':'-','color':r.color,'lw':.55,'alpha':.6})
a.set(xlim=(0,480),ylim=(.33,.70),xlabel='Mean geometric clusters at 0.5 Å',
      ylabel='MAT-R (Å; lower is better)')
a.set_xticks([100,200,300,400]);a.set_yticks([.35,.40,.45,.50,.55,.60,.65,.70])
offsets = {
    'core': {
        'ChEMBL3D-PB':(10,-20), 'RDKit raw':(10,3), 'RDKit minimized':(12,4),
        'Torsion raw':(-12,8), 'Torsion minimized':(10,8), 'LoQI':(-10,13),
        'Torsional Diffusion':(10,7), 'MCF drugs-L':(12,9),
        'NExT-Mol DMT-L':(12,-9), 'FlowR':(-12,-4), 'Qwen 1.7B FSQ':(-12,-17)},
    'ref': {
        'ChEMBL3D-PB':(10,-20), 'RDKit raw':(10,-28), 'RDKit minimized':(10,-9),
        'Torsion raw':(10,7), 'Torsion minimized':(10,6), 'LoQI':(10,-16),
        'Torsional Diffusion':(10,-20), 'MCF drugs-L':(10,12),
        'NExT-Mol DMT-L':(-10,-13), 'FlowR':(12,8), 'Qwen 1.7B FSQ':(10,-17)},
}
casf_styles={name:(color,marker) for name,color,marker in styles.values()}
annotations=[]
for ax, cohort, title, subtitle in zip(axes[1:], ['core','ref'],
        ['B   CASF core · 94 entries','C   CASF ref · 1,044 entries'],
        ['ChEMBL-count target','ChEMBL-count target']):
    ax.set_title(title,loc='left',pad=29,fontsize=11)
    ax.text(0,1.035,subtitle,transform=ax.transAxes,ha='left',fontsize=8,color='#626C78')
    ax.grid(True,color='#E7EBEE',linewidth=.6,zorder=0);ax.set_axisbelow(True)
    ax.tick_params(length=3,width=.6)
    d=points.loc[points.cohort.eq(cohort)]
    for r in d.itertuples():
        color,marker=casf_styles[r.method]
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

fig.text(.52,.035,'Axes differ across datasets. Dashed guides: stored ChEMBL3D-PB within each CASF panel.',ha='center',fontsize=8,color='#626C78')
def save(fig,name):
 for ext in ['png','pdf','svg']:
  fig.savefig(OUT/f'{name}.{ext}',dpi=300)
 svg=OUT/f'{name}.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
 plt.close(fig)
save(fig,'figure-diversity-proximity')
fig,axes=plt.subplots(1,2,figsize=(10,4.9))
fig.subplots_adjust(left=.09,right=.975,top=.73,bottom=.30,wspace=.30)
fig.suptitle('PLINDER: reference coverage and sample precision',x=.09,y=.975,ha='left',fontsize=12)
fig.text(.09,.915,'22 molecules · 1,067 cleaned references · supplied ensembles',fontsize=9,color='#626C78')
for ax in axes:
 ax.grid(color='#E7EBEE',lw=.6);ax.set_axisbelow(True);ax.tick_params(length=3,width=.6)
for ax,col,title in zip(axes,['recall','precision'],
                        ['A   Coverage of experimental conformations','B   Precision relative to observed conformations']):
    for method in curve_methods:
        name,color,marker=styles[methods[method]]
        d=curves.loc[curves.method.eq(method)]
        ax.plot(d.threshold,100*d[col],color=color,lw=2,
                ls=(0,(4,2)) if method == 'rdkit_random_raw_druglike' else '-',
                label=name)
        value=100*d.loc[d.threshold.eq(.75),col].iloc[0]
        ax.scatter([.75],[value],s=20,color=color,zorder=4,edgecolor='white',linewidth=.4)
    ax.axvline(.75,color='#626C78',ls=(0,(4,3)),lw=1,zorder=1)
    ax.text(.75,1.025,'0.75 Å',transform=ax.get_xaxis_transform(),ha='center',fontsize=9,color='#626C78')
    ax.set_title(title,loc='left',pad=27,fontsize=10)
    ax.set(xlim=(0,1.5),ylim=(0,100),xlabel='RMSD threshold (Å)')
    ax.set_xticks([0,.25,.5,.75,1,1.25,1.5],labels=['0','0.25','0.50','0.75','1.00','1.25','1.50'])
    ax.set_yticks([0,20,40,60,80,100])
axes[0].set_ylabel('Reference coverage, COV-R (%)')
axes[1].set_ylabel('Sample precision, COV-P (%)')
handles,labels=axes[0].get_legend_handles_labels()
fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.53,.085),ncol=4,
           frameon=False,fontsize=8.5,columnspacing=2.3,handlelength=2.7)
fig.text(.53,.028,'Dashed curves: RDKit raw.',
         ha='center',fontsize=7.7,color='#626C78')
save(fig,'supporting/plinder-thresholds')
sources=[PL/'summary.csv',PL/'threshold-curves.csv',DATA/'table1.csv',DATA/'table2.csv',PAPER/'sampling-data/palette.json',Path(__file__)]
outputs=[OUT/f'{name}.{ext}' for name in ['figure-diversity-proximity','supporting/plinder-thresholds'] for ext in ['png','pdf','svg']]
(HERE/'provenance.json').write_text(json.dumps({
 'inputs':{str(f.relative_to(PAPER)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources},
 'outputs':{str(f.relative_to(PAPER)):hashlib.sha256(f.read_bytes()).hexdigest() for f in outputs},
 'description':'Main Figure 1 combines PLINDER diversity/MAT-R and CASF diversity/Best RMSD. Supporting Figure S8 shows PLINDER COV-R/P versus threshold. All values unchanged.'
},indent=2)+'\n')
print('Rendered main Figure 1 and Supporting Figure S8 from unchanged source data.')
