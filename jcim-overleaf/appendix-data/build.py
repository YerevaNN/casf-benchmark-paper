"""Render the selected supporting analyses from current manuscript archives.

Run from any directory with the benchmark Python environment. No generation,
optimization, filtering, alignment, or manuscript compilation is performed.
"""
from pathlib import Path
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
OUT = PAPER / 'figures/supporting'
OUT.mkdir(exist_ok=True)
SOURCES = [
    'figure-data/candidates/core-entries.csv',
    'results-data/ref-panel/per-entry.csv',
    'results-data/three-dataset-tables/table3.csv',
    'results-data/size-flexibility/figure5-strata.csv',
    'plinder22-data/per-molecule.csv',
    'plinder22-data/qwen-dedup/per-molecule.csv',
    'results-data/druglike_means.csv',
    'results-data/qwen-size-representation/plinder-per-molecule.csv',
    'results-data/qwen-size-representation/08-main-size-and-representation-values.csv',
    'sampling-data/palette.json', 'plinder22-data/methods.json',
]
styles = json.loads((PAPER / 'sampling-data/palette.json').read_text())
keys = [k for k in styles if k != 'chembl3d_gt_pb']
REF = 'chembl3d_gt_pb'
QWEN = 'qwen_1p7b_v2_coordtoken_step41054'
OLD = 'qwen_0p6b_casf_dedup_binned_pm12_step64000'
INK = '#29323C'
plt.rcParams.update({
    'font.family':'DejaVu Sans', 'font.size':9, 'axes.titlesize':11,
    'axes.labelsize':9, 'xtick.labelsize':8, 'ytick.labelsize':8,
    'axes.edgecolor':'#8B929A', 'axes.labelcolor':INK, 'text.color':INK,
    'xtick.color':'#535D68', 'ytick.color':'#535D68', 'axes.linewidth':.65,
    'axes.spines.top':False, 'axes.spines.right':False,
    'pdf.fonttype':42, 'svg.fonttype':'none', 'svg.hashsalt':'curated-supporting-20261105',
    'savefig.facecolor':'white', 'figure.facecolor':'white',
})
outputs = []

def save(fig, stem):
    for extension in ['pdf','png','svg']:
        path = OUT / f'{stem}.{extension}'
        metadata = {'CreationDate':None,'ModDate':None} if extension=='pdf' else {'Date':None} if extension=='svg' else None
        fig.savefig(path, dpi=250, metadata=metadata)
        if extension=='svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
        outputs.append(path)
    plt.close(fig)

# Pair individual hits before summarizing, retaining the complete panel denominator.
table = pd.read_csv(PAPER / SOURCES[2])
strata = pd.read_csv(PAPER / SOURCES[3])
casf = {c:pd.read_csv(PAPER / source) for c,source in zip(['core','ref'],SOURCES[:2])}
changes, composition = [], []
for cohort, n in [('core',94),('ref',1044)]:
    d = casf[cohort]
    base = d[d.method.eq(REF)].set_index('mol_id')
    assert base.index.is_unique and len(base)==n
    base_hit = base.casf_hit_0p75.fillna(0)
    for key in keys:
        part = d[d.method.eq(key+'_fixed')].set_index('mol_id')
        assert part.index.is_unique and set(part.index)==set(base.index)
        hit = part.reindex(base.index).casf_hit_0p75.fillna(0)
        assert hit.isin([0,1]).all() and base_hit.isin([0,1]).all()
        gained=int(((hit==1)&(base_hit==0)).sum())
        lost=int(((hit==0)&(base_hit==1)).sum())
        shared=int(((hit==1)&(base_hit==1)).sum())
        expected=table[(table.cohort==cohort)&(table.method==styles[key][0])&(table.budget=='fixed')].iloc[0]
        np.testing.assert_allclose(hit.mean()*100, expected.recovery, atol=1e-12)
        assert gained-lost==int(hit.sum()-base_hit.sum())
        changes.append(dict(cohort=cohort,method=key,name=styles[key][0],n=n,gained=gained,lost=lost,shared_hits=shared,shared_misses=int(((hit==0)&(base_hit==0)).sum()),net=gained-lost))
    size = pd.cut(base.heavy_atoms,[0,20,30,40,np.inf],right=False,labels=['<20','20-29','30-39','40+'])
    rot = pd.cut(base.rotatable_bonds,[-1,3,6,8,np.inf],labels=['0-3','4-6','7-8','9+'])
    assert size.notna().all() and rot.notna().all()
    counts=pd.crosstab(rot,size,dropna=False)
    assert counts.to_numpy().sum()==n
    for descriptor,series in [('size',counts.sum(axis=0)),('rotors',counts.sum(axis=1))]:
        expected=strata[(strata.cohort==cohort)&(strata.key==REF)&(strata.descriptor==descriptor)].set_index('group')
        np.testing.assert_array_equal(series.to_numpy(),expected.loc[series.index,'n'])
    for rotor in counts.index:
        for heavy in counts.columns:
            composition.append(dict(cohort=cohort,size=str(heavy),rotors=str(rotor),count=int(counts.loc[rotor,heavy]),percent=100*counts.loc[rotor,heavy]/n))
changes=pd.DataFrame(changes); changes.to_csv(HERE/'recovery-gains-losses.csv',index=False)
composition=pd.DataFrame(composition);composition.to_csv(HERE/'cohort-composition.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(10,5.6),sharey=True)
fig.subplots_adjust(left=.20,right=.98,bottom=.17,top=.87,wspace=.26)
fig.suptitle('Recovery gained and lost relative to ChEMBL3D',x=.20,ha='left',y=.98,fontsize=13)
order=keys # shared established method order, Qwen last
for ax,cohort,letter in zip(axes,['core','ref'],['A','B']):
    d=changes[changes.cohort==cohort].set_index('method').loc[order]
    maximum=max(d.gained.max(),d.lost.max()); pad=max(1,.035*maximum)
    for y,(key,r) in enumerate(d.iterrows()):
        color=styles[key][1]
        ax.barh(y,-r.lost,color='white',edgecolor=color,hatch='///',height=.63,lw=.9)
        ax.barh(y,r.gained,color=color,height=.63)
        ax.text(-r.lost-pad,y,str(r.lost),ha='right',va='center',fontsize=8)
        ax.text(r.gained+pad,y,str(r.gained),ha='left',va='center',fontsize=8)
    ax.set_title(f'{letter}   CASF {cohort} ({int(d.n.iloc[0]):,})',loc='left',pad=12)
    ax.axvline(0,color=INK,lw=.8)
    ax.set_xlim(-d.lost.max()-5*pad,d.gained.max()+5*pad)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{abs(x):g}'))
    ax.set_xlabel('Stored only  ←  Ligands  →  Generated only',labelpad=10)
    ax.grid(axis='x',color='#E7EBEE',lw=.6);ax.set_axisbelow(True)
axes[0].set_yticks(range(len(order)),[styles[k][0] for k in order]);axes[0].invert_yaxis()
save(fig,'recovery-gains-losses')

fig,axes=plt.subplots(1,2,figsize=(8.8,4.3))
fig.subplots_adjust(left=.10,right=.95,bottom=.20,top=.78,wspace=.42)
fig.suptitle('Joint distribution of ligand size and flexibility',x=.10,ha='left',y=.98,fontsize=12)
size_order=['<20','20-29','30-39','40+'];rot_order=['0-3','4-6','7-8','9+']
for ax,cohort,letter in zip(axes,['core','ref'],['A','B']):
    d=composition[composition.cohort==cohort]
    counts=d.pivot(index='rotors',columns='size',values='count').loc[rot_order,size_order]
    perc=d.pivot(index='rotors',columns='size',values='percent').loc[rot_order,size_order]
    im=ax.imshow(perc,cmap='Blues',vmin=0,vmax=45,aspect='auto')
    for y,x in np.ndindex(counts.shape):
        ax.text(x,y,f'{counts.iloc[y,x]}\n({perc.iloc[y,x]:.1f}%)',ha='center',va='center',fontsize=9,color='white' if perc.iloc[y,x]>27 else INK)
    ax.set_title(f'{letter}   CASF {cohort} ({counts.to_numpy().sum():,})',loc='left',pad=11)
    ax.set_xticks(range(4),['<20','20–29','30–39','≥40']);ax.set_yticks(range(4),['0–3','4–6','7–8','≥9'])
    ax.set_xlabel('Heavy atoms');ax.set_ylabel('Rotatable bonds')
    ax.tick_params(length=0)
    for spine in ax.spines.values():spine.set_visible(False)
save(fig,'cohort-composition')

# All eleven Table 4 methods, including the new independent torsion pools.
pl=pd.concat([pd.read_csv(PAPER/SOURCES[4]),pd.read_csv(PAPER/SOURCES[5])],ignore_index=True)
labels=json.loads((PAPER/'plinder22-data/methods.json').read_text())
labels[OLD]='Qwen 0.6B CASF-dedup'
means=pd.read_csv(PAPER/'results-data/druglike_means.csv').set_index('label')
methods=list(labels)
assert len(methods)==11 and len(pl)==242 and not pl.duplicated(['label','smiles']).any()
names=pl.drop_duplicates('smiles').set_index('smiles').name.sort_values()
assert len(names)==22 and names.is_unique
for method in methods:
    d=pl[pl.label==method].set_index('smiles').loc[names.index]
    assert d.num_true_confs.sum()==1067
    for field,target in [('cov_r_0p75','cov_r_075'),('cov_p_0p75','cov_p_075')]:
        np.testing.assert_allclose(d[field].mean(),means.loc[method,target],atol=1e-12)
    np.testing.assert_allclose(d.n_pb_pass.sum()/d.n_generated.sum(),means.loc[method,'pooled_pb_pass'],atol=1e-12)
pl[['label','smiles','name','num_true_confs','n_generated','n_pb_pass','cov_r_0p75','cov_p_0p75']].to_csv(HERE/'plinder-by-molecule.csv',index=False)
short_labels=['RDKit\nraw','RDKit\nminimized','Torsion\nraw','Torsion\nminimized','LoQI','Torsional\nDiffusion','MCF\ndrugs-L','NExT-Mol\nDMT-L','FlowR','Qwen 1.7B\nFSQ','Qwen 0.6B\nCASF-dedup']

def heatmap(ax,values,column_names,show_names=True):
    im=ax.imshow(values,cmap='Blues',norm=Normalize(0,100),aspect='auto')
    for y,x in np.ndindex(values.shape):
        value=values[y,x]
        ax.text(x,y,f'{value:.0f}',ha='center',va='center',fontsize=7.8,color='white' if value>=60 else INK)
    ax.set_xticks(range(len(column_names)),column_names,rotation=55,ha='right',fontsize=8)
    ax.set_yticks(range(22),names.tolist() if show_names else [],fontsize=9)
    ax.tick_params(which='both',length=0)
    ax.set_xticks(np.arange(-.5,len(column_names)),minor=True)
    ax.set_yticks(np.arange(-.5,22),minor=True)
    ax.grid(which='minor',color='white',lw=.55)
    for spine in ax.spines.values():spine.set_visible(False)
    return im

fig,axes=plt.subplots(1,2,figsize=(12.4,8.2))
fig.subplots_adjust(left=.16,right=.985,bottom=.21,top=.88,wspace=.11)
fig.suptitle('PLINDER: coverage and precision for each molecule',x=.16,ha='left',y=.98,fontsize=13)
fig.text(.16,.945,'22 molecules · 1,067 references · RMSD ≤ 0.75 Å',color='#626C78',fontsize=10)
for i,(ax,metric,title) in enumerate(zip(axes,['cov_r_0p75','cov_p_0p75'],['A   Reference coverage (COV-R)','B   Sample precision (COV-P)'])):
    values=pl.pivot(index='smiles',columns='label',values=metric).loc[names.index,methods].to_numpy()*100
    im=heatmap(ax,values,short_labels,i==0);ax.set_title(title,loc='left',pad=12)
cax=fig.add_axes([.43,.045,.32,.015]);fig.colorbar(im,cax=cax,orientation='horizontal',label='Percentage (%)')
save(fig,'plinder-coverage-precision')

# Alongside Table 4 validity, show all four v2 Qwen configurations on the same molecules.
v2=pd.read_csv(PAPER/'results-data/qwen-size-representation/plinder-per-molecule.csv')
# Column names are verified against the archived four-configuration export.
qorder=['binned06','binned17','fsq06','fsq17']
qkey='model'
assert set(v2[qkey])==set(qorder)
qmeans=pd.read_csv(PAPER/'results-data/qwen-size-representation/08-main-size-and-representation-values.csv')
for key in qorder:
    d=v2[v2[qkey]==key]
    assert len(d)==22 and set(d.smiles)==set(names.index)
    exp=qmeans[(qmeans.panel=='C')&(qmeans.representation==('FSQ' if key.startswith('fsq') else 'Binned'))&(qmeans['size']==(0.6 if key.endswith('06') else 1.7))].value.iloc[0]
    np.testing.assert_allclose(100*d.n_pb_pass.sum()/d.n_generated.sum(),exp,atol=1e-10)
v2[['model','label','smiles','name','n_generated','n_pb_pass','pb_pass_rate']].to_csv(HERE/'qwen-validity-source.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(12.4,8.2),gridspec_kw={'width_ratios':[11,4]})
fig.subplots_adjust(left=.16,right=.985,bottom=.21,top=.88,wspace=.12)
fig.suptitle('PLINDER: structural validity for each molecule',x=.16,ha='left',y=.98,fontsize=13)
fig.text(.16,.945,'PoseBusters pass percentages · validity assessed separately from matching',color='#626C78',fontsize=10)
values=pl.pivot(index='smiles',columns='label',values='pb_pass_rate').loc[names.index,methods].to_numpy()*100
im=heatmap(axes[0],values,short_labels);axes[0].set_title('A   Generators in Table 4',loc='left',pad=12)
values=v2.pivot(index='smiles',columns=qkey,values='pb_pass_rate').loc[names.index,qorder].to_numpy()*100
heatmap(axes[1],values,['0.6B\nbinned','1.7B\nbinned','0.6B\nFSQ','1.7B\nFSQ'],False)
axes[1].set_title('B   Four Qwen models',loc='left',pad=12)
cax=fig.add_axes([.43,.045,.32,.015]);fig.colorbar(im,cax=cax,orientation='horizontal',label='PoseBusters pass (%)')
save(fig,'plinder-validity')

# Retained ensemble sizes are supporting evidence not shown in the main budget table.
lines=[r'% Rebuild with appendix-data/build.py.',r'\begin{table}[htbp]',r'\centering\small',r'\caption{Mean retained ensemble sizes at the two CASF candidate targets.}',r'\label{tab:supporting-ensemble-counts}',r'\begin{tabular}{lrrrr}',r'\toprule',r' & \multicolumn{2}{c}{CASF core (94)} & \multicolumn{2}{c}{CASF ref (1,044)} \\',r'\cmidrule(lr){2-3}\cmidrule(l){4-5}',r'Method & ChEMBL count & 1,000 & ChEMBL count & 1,000 \\',r'\midrule']
for key in keys+[REF]:
    label=styles[key][0];cells=[]
    if key==REF:lines.append(r'\hdashline')
    for cohort in ['core','ref']:
        for budget in ['chembl_count','fixed']:
            if key==REF and budget=='fixed':cells.append('---');continue
            row=table[(table.method==label)&(table.cohort==cohort)&(table.budget==budget)].iloc[0]
            cells.append(f'{row.mean_retained:.1f}')
    lines.append(label+' & '+' & '.join(cells)+r' \\')
lines += [r'\bottomrule',r'\end{tabular}',r'\par\smallskip',r'\begin{minipage}{\linewidth}',r'\footnotesize Counts are measured after PoseBusters filtering and averaged over all entries, including empty outputs. ChEMBL3D-PB is the same stored ensemble at both targets and is shown once. Rejected conformers are not replaced.',r'\end{minipage}',r'\end{table}']
target=PAPER/'tables/supporting-ensemble-counts.tex';target.write_text('\n'.join(lines)+'\n');outputs.append(target)
for name in ['recovery-gains-losses.csv','cohort-composition.csv','plinder-by-molecule.csv','qwen-validity-source.csv']:
    outputs.append(HERE/name)
(HERE/'provenance.json').write_text(json.dumps({
    'scope':'Selected supporting figures S3–S6 and retained-count Table S1. Current core94/ref1044/PLINDER22; primary Qwen v2 step41054; independent torsion baselines.',
    'inputs':{name:hashlib.sha256((PAPER/name).read_bytes()).hexdigest() for name in SOURCES},
    'outputs':{str(path.relative_to(PAPER)):hashlib.sha256(path.read_bytes()).hexdigest() for path in outputs},
    'checks':['All CASF hit means reproduce Table 3','Gained minus lost equals the net recovery difference','Joint descriptor margins reproduce Figure 5','All eleven PLINDER methods reproduce Table 4','Four v2 pooled validity values reproduce Figure 7'],
},indent=2)+'\n')
print('Built four supporting figures and the retained-count table from current records.')
