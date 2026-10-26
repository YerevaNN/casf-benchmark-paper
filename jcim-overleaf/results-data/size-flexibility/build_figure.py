"""Rebuild Figure 5 from archived CASF core and ref strata."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path(__file__).resolve().parent
PAPER = OUT.parents[1]
styles = json.loads((PAPER/'sampling-data/palette.json').read_text())
name_to_key = {v[0]:k for k,v in styles.items()}
ink = '#29323C'
TITLES = {'core':'CASF core (94)', 'ref':'CASF ref (1,236)'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.4,'pdf.fonttype':42,
                    'svg.fonttype':'none','svg.hashsalt':'size-flexibility-three-datasets',
                    'text.color':ink,'axes.labelcolor':ink,'savefig.facecolor':'white'})

# Extend the existing size/flexibility analysis from the same archived CASF exports.
order=pd.read_csv(PAPER/'results-data/three-dataset-tables/table1.csv').method.drop_duplicates().tolist()
keys=[name_to_key[m] for m in order]
size=pd.read_csv(OUT/'casf-size-source.csv');flex=pd.read_csv(OUT/'casf-flexibility-source.csv')
strata=[]
groups={'size':['<20','20-29','30-39','40+'],'rotors':['0-3','4-6','7-8','9+']}
for cohort,total in [('core',94),('ref',1236)]:
    for key in keys:
        method=key if key=='chembl3d_gt_pb' else key+'_fixed'
        for descriptor,source,col in [('size',size,'heavy_atoms'),('rotors',flex,'rotors')]:
            rows=source.loc[source['set'].eq(cohort)&source.method.eq(method)]
            assert len(rows)==4 and rows.n.sum()==total,(cohort,method,descriptor)
            for r in rows.to_dict('records'):
                hits=r.get('hits',r['hit_075']*r['n']);assert np.isclose(hits,round(hits))
                strata.append(dict(cohort=cohort,key=key,descriptor=descriptor,group=r[col],n=int(r['n']),hits=int(round(hits)),value=100*r['hit_075'],metric='Hit@0.75'))
casf=pd.DataFrame(strata)
old=pd.read_csv(PAPER/'results-data/strata.csv')
check=casf.query("cohort == 'core'").merge(old,on=['key','descriptor','group'])
assert len(check)==len(old)==88
np.testing.assert_allclose(check.value,check.hit)
np.testing.assert_array_equal(check.n_x,check.n_y)
np.testing.assert_array_equal(check.hits_x,check.hits_y)

all_strata=casf
all_strata.to_csv(OUT/'figure5-strata.csv',index=False)

def sizefigure(cohorts):
    fig,axes=plt.subplots(len(cohorts),2,figsize=(7.2,5.8),squeeze=False)
    fig.subplots_adjust(left=.24,right=.98,top=.91,bottom=.17,hspace=.64,wspace=.23)
    labels={'size':['<20','20–29','30–39','≥40'],'rotors':['0–3','4–6','7–8','≥9']}
    for i,cohort in enumerate(cohorts):
        for j,descriptor in enumerate(groups):
            ax=axes[i,j];d=all_strata.loc[all_strata.cohort.eq(cohort)&all_strata.descriptor.eq(descriptor)].set_index(['key','group'])
            matrix=np.full((len(keys),4),np.nan)
            counts=[]
            for g in groups[descriptor]:counts.append(int(d.xs(g,level='group').n.iloc[0]))
            for y,key in enumerate(keys):
                for x,g in enumerate(groups[descriptor]):
                    if (key,g) in d.index:matrix[y,x]=d.loc[(key,g),'value']
            cmap=plt.get_cmap('Blues').copy();cmap.set_bad('#EEF0F2')
            mesh=ax.imshow(matrix,cmap=cmap,vmin=0,vmax=100,aspect='auto')
            for y,key in enumerate(keys):
                for x,g in enumerate(groups[descriptor]):
                    value=matrix[y,x]
                    if np.isnan(value):text='—'
                    else:
                        r=d.loc[(key,g)];text=f'{int(r.hits)}/{int(r.n)}'
                    ax.text(x,y,text,ha='center',va='center',fontsize=7.2,color='white' if value>=65 else ink)
            ax.set_xticks(range(4),[f'{label}{"*" if 0<n<5 else ""}\n(n={n})' for label,n in zip(labels[descriptor],counts)],fontsize=7.4)
            ax.set_yticks(range(len(keys)),[styles[k][0] for k in keys] if j==0 else [],fontsize=7.4)
            ax.tick_params(which='both',length=0)
            ax.set_xticks(np.arange(-.5,4),minor=True);ax.set_yticks(np.arange(-.5,len(keys)),minor=True);ax.grid(which='minor',color='white',lw=.9)
            ax.set_title(f'{TITLES[cohort]} · Hit@0.75\n'+('Heavy atoms' if j==0 else 'Rotatable bonds'),loc='left',fontsize=8.2,pad=9)
            for x,n in enumerate(counts):
                if 0<n<5:ax.add_patch(Rectangle((x-.5,-.5),1,len(keys),fill=False,edgecolor='#CE6428',ls='--',lw=1,zorder=5,clip_on=False))
            for spine in ax.spines.values():spine.set_visible(False)
    cax=fig.add_axes([.40,.065,.42,.012]);fig.colorbar(mesh,cax=cax,orientation='horizontal',label='Recovery (%)')
    return fig

fig = sizefigure(['core','ref'])
for extension in ['pdf','svg','png']:
    metadata = {'CreationDate':None, 'ModDate':None} if extension=='pdf' else {'Date':None} if extension=='svg' else None
    target = PAPER/'figures'/f'figure-4-size-flexibility.{extension}'
    fig.savefig(target,dpi=300,metadata=metadata)
    if extension == 'svg':
        target.write_text('\n'.join(line.rstrip() for line in target.read_text().splitlines())+'\n')
plt.close(fig)
sources = [OUT/'casf-size-source.csv', OUT/'casf-flexibility-source.csv',
           PAPER/'results-data/strata.csv',
           PAPER/'sampling-data/palette.json', PAPER/'results-data/three-dataset-tables/table1.csv']
(OUT/'provenance.json').write_text(json.dumps({
    'sources':{str(p.relative_to(PAPER)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    'datasets':['CASF core94','CASF ref1236'],
    'aggregation':'CASF: Hit@0.75 at PB-filtered fixed target; stored baseline unchanged.',
    'qwen':'Original Qwen 1.7B FSQ step47023; no retrained outputs substituted.',
    'checks':['All 88 existing CASF core cells reproduced','CASF group denominators sum to 94/1236']
},indent=2)+'\n')
print('Built Figure 5 for CASF core/ref; all 88 original core cells reproduced.')
