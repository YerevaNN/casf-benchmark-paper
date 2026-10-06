"""Preview threshold-dependent PLINDER coverage from existing distance caches."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path('/mnt/weka/mbedrosian/code/casf-benchmark')
OUT=Path(__file__).resolve().parent
CACHE=ROOT/'data/results/cache/druglike_rmsd_matrices'
STYLE=json.loads((ROOT/'docs/sampling_analysis/palette.json').read_text())
METHODS={
'qwen_1p7b_fsq_bigdata_step47023':'qwen_1p7b_fsq_bigdata_step47023',
'loqi_druglike':'loqi_raw',
'torsional_diffusion_druglike':'torsional_diffusion_raw',
'mcf_drugs_l_druglike':'mcf_drugs_l_raw',
'flowr_druglike':'flowr_raw',
'nextmol_dmt_l_druglike':'nextmol_dmt_l_raw',
}
thresholds=np.linspace(0,1.5,301)
assert thresholds[150]==.75
rows=[];audit=[];hashes={};identities=None
archived=pd.read_csv(ROOT/'docs/results_tables/druglike_selected.csv')
table=pd.read_csv(ROOT/'docs/results_tables/druglike_means.csv').set_index('label')
for method,key in METHODS.items():
    path=CACHE/(method+'.npz');hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    recalls=[];precisions=[];old=archived[archived.label.eq(method)].set_index('smiles')
    with np.load(path,allow_pickle=True) as cache:
        smiles=cache['smiles'].tolist();assert len(smiles)==23
        if identities is None:identities=set(smiles)
        assert set(smiles)==identities==set(old.index)
        for i,smi in enumerate(smiles):
            matrix=cache[f'm{i}']
            expected=old.loc[smi]
            assert matrix.shape==(int(expected.num_true_confs),int(expected.num_gen_confs)),(method,smi,matrix.shape)
            finite=np.isfinite(matrix);assert np.all(matrix[finite]>=0)
            matrix=np.where(finite,matrix,np.inf)
            reference_min=matrix.min(axis=1);generated_min=matrix.min(axis=0)
            # No all-failed rows/columns: all original records stay in the denominators.
            assert np.isfinite(reference_min).all() and np.isfinite(generated_min).all()
            r=(reference_min[:,None]<thresholds).mean(axis=0)
            p=(generated_min[:,None]<thresholds).mean(axis=0)
            recalls.append(r);precisions.append(p)
            audit.append(dict(method=method,smiles=smi,n_reference=len(reference_min),n_generated=len(generated_min),
              cached_recall_075=r[150],table_recall_075=expected.cov_r_075,
              cached_precision_075=p[150],table_precision_075=expected.cov_p_075))
    for t,r,p in zip(thresholds,np.mean(recalls,axis=0),np.mean(precisions,axis=0)):
        rows.append(dict(method=method,threshold=t,recall=r,precision=p))
curves=pd.DataFrame(rows)
curves.to_csv(OUT/'threshold_curves.csv',index=False)
pd.DataFrame(audit).to_csv(OUT/'cache_table_comparison.csv',index=False)
for method in METHODS:
    c=curves[curves.method.eq(method)]
    assert (np.diff(c.recall)>=0).all() and (np.diff(c.precision)>=0).all()
    assert c[['recall','precision']].ge(0).all().all() and c[['recall','precision']].le(1).all().all()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
 'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,
 'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C',
 'xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
 'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})
fig,axes=plt.subplots(1,2,figsize=(10.8,5),sharex=True,sharey=True)
fig.subplots_adjust(left=.075,right=.975,top=.77,bottom=.28,wspace=.24)
for ax,col,title in zip(axes,['recall','precision'],['A   Coverage of experimental conformations','B   Precision relative to observed conformations']):
    for method,key in METHODS.items():
        label,color,_=STYLE[key];c=curves[curves.method.eq(method)]
        ax.plot(c.threshold,c[col]*100,label=label,color=color,lw=2)
        value=c.loc[c.threshold.eq(.75),col].iloc[0]*100
        ax.scatter([.75],[value],s=20,color=color,zorder=4,edgecolor='white',linewidth=.4)
    ax.axvline(.75,color='#626C78',ls=(0,(4,3)),lw=1,zorder=1)
    ax.text(.75,1.025,'0.75 Å',transform=ax.get_xaxis_transform(),ha='center',fontsize=9,color='#626C78')
    ax.set_title(title,loc='left',pad=27)
    ax.set_xlim(0,1.5);ax.set_ylim(0,100)
    ax.set_xticks([0,.25,.5,.75,1,1.25,1.5],labels=['0','0.25','0.50','0.75','1.00','1.25','1.50'])
    ax.set_yticks([0,20,40,60,80,100]);ax.grid(color='#E7EBEE',lw=.6);ax.set_axisbelow(True)
    ax.set_xlabel('RMSD threshold (Å)')
axes[0].set_ylabel('Reference coverage, COV-R (%)')
axes[1].set_ylabel('Sample precision, COV-P (%)')
fig.suptitle('Coverage and precision across matching tolerances',x=.075,y=.975,ha='left',fontsize=13)
fig.text(.075,.913,'PLINDER: 23 molecules · equal weight per molecule · supplied conformer ensembles',fontsize=9,color='#626C78')
handles,labels=axes[0].get_legend_handles_labels()
fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.52,.077),ncol=3,frameon=False,fontsize=9,columnspacing=2,handlelength=2.7)
fig.text(.5,.02,'Preview from cached distances; strict RMSD < threshold. Common validity filtering remains pending.',ha='center',fontsize=8,color='#626C78')
for ext in ['png','pdf','svg']:fig.savefig(OUT/('plinder-threshold-preview.'+ext),dpi=250)
plt.close(fig)
(OUT/'provenance.json').write_text(json.dumps({'cache_sha256':hashes,'n_molecules':23,'threshold_rule':'strict <','threshold_range':[0,1.5],'threshold_step':.005,'aggregation':'equal-weight molecule mean','validity':'supplied pools; no common PB filtering','alignment_recomputed':False,'table_unchanged':True},indent=2)+'\n')
print(curves[curves.threshold.isin([.25,.5,.75,1,1.5])].round(4).to_string(index=False))
print('Validated all six caches: 23 matching molecular identities and matching reference/generated counts; finite minima; monotonic bounded curves.')
