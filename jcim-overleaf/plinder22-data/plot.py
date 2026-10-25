"""Render the combined PLINDER figure from the validated release archive."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
OUT = PAPER / 'figures'
styles = json.loads((PAPER / 'sampling-data/palette.json').read_text())
methods = {
    'qwen_1p7b_fsq_bigdata_step47023':'qwen_1p7b_fsq_bigdata_step47023',
    'loqi_druglike':'loqi_raw', 'torsional_diffusion_druglike':'torsional_diffusion_raw',
    'mcf_drugs_l_druglike':'mcf_drugs_l_raw', 'flowr_druglike':'flowr_raw',
    'nextmol_dmt_l_druglike':'nextmol_dmt_l_raw',
    'rdkit_random_raw_druglike':'rdkit_random_raw', 'torsion_raw_druglike':'torsion_raw',
    'rdkit_random_minimized_22_druglike':'rdkit_random_minimized',
    'torsion_minimized_22_druglike':'torsion_minimized',
}
classical = {'rdkit_random_raw_druglike','torsion_raw_druglike'}
minimized = {'rdkit_random_minimized_22_druglike','torsion_minimized_22_druglike'}
seeded = {'torsion_raw_druglike','torsion_minimized_22_druglike'}
summary = pd.read_csv(HERE / 'summary.csv').set_index('label')
points = pd.DataFrame([dict(method=m, name=styles[k][0], color=styles[k][1],
    marker=styles[k][2], clusters_0p5=summary.loc[m, 'mean_clusters_0p5'],
    mat_r=summary.loc[m, 'mat_r_mean']) for m,k in methods.items()])
points.to_csv(HERE / 'scatter-points.csv', index=False)
curves = pd.read_csv(HERE / 'threshold-curves.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,
    'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
    'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C',
    'xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,
    'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})
fig,axes=plt.subplots(1,3,figsize=(17.5,6),gridspec_kw={'width_ratios':[1.12,1,1]})
fig.subplots_adjust(left=.055,right=.985,top=.76,bottom=.30,wspace=.28)
fig.suptitle('PLINDER: geometric diversity, reference coverage and sample precision',
             x=.055,y=.975,ha='left',fontsize=13)
fig.text(.055,.923,'22 molecules · 1,067 cleaned references · multiple observed bound conformations · equal weight per molecule · supplied ensembles',
         fontsize=9,color='#626C78')
for ax in axes:
    ax.grid(color='#E7EBEE',lw=.6);ax.set_axisbelow(True)
    ax.tick_params(length=3,width=.6)
a=axes[0]
a.set_title('A   Diversity versus proximity',loc='left',pad=27)
offsets={'RDKit minimized':(-8,7),'Torsion minimized':(10,-12),'RDKit raw':(10,7),'Torsion raw':(-10,12),'LoQI':(10,8),
         'Torsional Diffusion':(10,-15),'MCF drugs-L':(10,17),
         'NExT-Mol DMT-L':(10,5),'FlowR':(-12,7),'Qwen 1.7B FSQ':(-12,-19)}
for r in points.itertuples():
    a.scatter(r.clusters_0p5,r.mat_r,s=65,marker=r.marker,color=r.color,
              edgecolor='white',linewidth=.65,zorder=5)
    dx,dy=offsets[r.name]
    a.annotate(r.name+('†' if r.method in seeded else ''),
               (r.clusters_0p5,r.mat_r),xytext=(dx,dy),textcoords='offset points',
               ha='right' if dx<0 else 'left',va='bottom',fontsize=8,color=r.color,
               arrowprops={'arrowstyle':'-','color':r.color,'lw':.55,'alpha':.6})
a.set(xlim=(0,480),ylim=(.33,.60),xlabel='Mean geometric clusters at 0.5 Å',
      ylabel='MAT-R (Å; lower is better)')
a.set_xticks([100,200,300,400]);a.set_yticks([.35,.40,.45,.50,.55,.60])
for ax,col,title in zip(axes[1:],['recall','precision'],
                        ['B   Coverage of experimental conformations','C   Precision relative to observed conformations']):
    for method,key in methods.items():
        name,color,marker=styles[key]
        d=curves.loc[curves.method.eq(method)]
        ax.plot(d.threshold,100*d[col],color=color,lw=2,
                ls=(0,(1,1.5)) if method in minimized else (0,(4,2)) if method in classical else '-',
                label=name+('†' if method in seeded else ''))
        value=100*d.loc[d.threshold.eq(.75),col].iloc[0]
        ax.scatter([.75],[value],s=20,color=color,zorder=4,edgecolor='white',linewidth=.4)
    ax.axvline(.75,color='#626C78',ls=(0,(4,3)),lw=1,zorder=1)
    ax.text(.75,1.025,'0.75 Å',transform=ax.get_xaxis_transform(),ha='center',fontsize=9,color='#626C78')
    ax.set_title(title,loc='left',pad=27,fontsize=10)
    ax.set(xlim=(0,1.5),ylim=(0,100),xlabel='RMSD threshold (Å)')
    ax.set_xticks([0,.25,.5,.75,1,1.25,1.5],labels=['0','0.25','0.50','0.75','1.00','1.25','1.50'])
    ax.set_yticks([0,20,40,60,80,100])
axes[1].set_ylabel('Reference coverage, COV-R (%)')
axes[2].set_ylabel('Sample precision, COV-P (%)')
handles,labels=axes[1].get_legend_handles_labels()
fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.53,.135),ncol=5,
           frameon=False,fontsize=8.5,columnspacing=2.3,handlelength=2.7)
fig.text(.52,.093,'Full supplied pools for clustering, MAT-R and curves. RMSD ≤ threshold; no PB filter. Equal weight per molecule.',
         ha='center',fontsize=7.7,color='#626C78')
fig.text(.52,.052,'† Both torsion baselines start from an experimental reference. Classical curves: dashed, raw; dotted, minimized.',
         ha='center',fontsize=7.7,color='#626C78')
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'figure-plinder-combined.{ext}',dpi=300)
plt.close(fig)

# Keep generated vector files clean in source control.
svg = PAPER / 'figures/figure-plinder-combined.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
