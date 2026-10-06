"""Build the working Results figures from archived, traceable measurements."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
TABLES=ROOT/'docs/results_tables'
Q='qwen_1p7b_fsq_bigdata_step47023'
STYLE={
 'chembl3d_gt_pb':('ChEMBL3D-PB','#333A43','h'),
 'rdkit_random_raw':('RDKit raw','#6F7782','o'),
 'rdkit_random_minimized':('RDKit minimized','#6F7782','s'),
 'torsion_raw':('Torsion raw','#987049','o'),
 'torsion_minimized':('Torsion minimized','#987049','s'),
 'loqi_raw':('LoQI','#CE6428','s'),
 'torsional_diffusion_raw':('Torsional Diffusion','#AA589B','D'),
 'mcf_drugs_l_raw':('MCF drugs-L','#998000','v'),
 'nextmol_dmt_l_raw':('NExT-Mol DMT-L','#399BBB','P'),
 'flowr_raw':('FlowR','#00866B','^'),
 Q:('Qwen 1.7B FSQ','#0072B2','o'),
}
KEYS=list(STYLE)
DKEY={'loqi_raw':'loqi_druglike','torsional_diffusion_raw':'torsional_diffusion_druglike',
 'mcf_drugs_l_raw':'mcf_drugs_l_druglike','nextmol_dmt_l_raw':'nextmol_dmt_l_druglike',
 'flowr_raw':'flowr_druglike',Q:Q}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.titleweight':'medium','axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C',
 'xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
 'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})
casf=pd.read_csv(TABLES/'casf_selected.csv')
drug=pd.read_csv(TABLES/'druglike_selected.csv')
means=pd.read_csv(TABLES/'druglike_means.csv').set_index('label')
inputs=[TABLES/x for x in ['casf_selected.csv','druglike_selected.csv','druglike_means.csv']]
figures=[]
def row(key,tier):
 method=key if key=='chembl3d_gt_pb' else key+'_'+tier
 r=casf[casf.method.eq(method)]
 assert len(r)==1
 return r.iloc[0]
def axis(ax,title,grid='both'):
 ax.set_title(title,loc='left',pad=13)
 ax.grid(True,axis=grid,color='#E7EBEE',linewidth=.6,zorder=0)
 ax.set_axisbelow(True)
 ax.tick_params(length=3,width=.6)
def point(ax,key,x,y,size=65):
 _,color,marker=STYLE[key]
 ax.scatter(x,y,s=size,marker=marker,color=color,edgecolor='white',linewidth=.65,zorder=5)
def save(fig,name,caption):
 fig.savefig(OUT/(name+'.pdf'))
 fig.savefig(OUT/(name+'.svg'))
 fig.savefig(OUT/(name+'.png'),dpi=350)
 figures.append((name,caption))
 plt.close(fig)
def footer(fig,text):
 fig.text(.5,.025,text,ha='center',fontsize=7.5,color='#626C78')

# 1. Cohort construction. No invented counts for individual exclusion steps.
fig,ax=plt.subplots(figsize=(7.2,4.65));fig.subplots_adjust(left=.025,right=.975,top=.86,bottom=.08)
ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
def box(x,y,w,h,title,detail,color='#0072B2'):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.01,rounding_size=0.015',
  edgecolor=color,facecolor='#F5F8FA',linewidth=1))
 ax.text(x+w/2,y+h*.69,title,ha='center',va='center',fontsize=10,fontweight='medium',color=color)
 ax.text(x+w/2,y+h*.30,detail,ha='center',va='center',fontsize=8,linespacing=1.4)
def arrow(a,b):
 ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'-|>','color':'#8B929A','lw':1.1,'shrinkA':3,'shrinkB':3})
box(.015,.74,.255,.19,'CASF-2016','285 crystal complexes')
box(.015,.46,.255,.19,'ChEMBL3D','Stored conformer ensembles')
box(.365,.60,.255,.23,'Core panel: 94','Identity + stereochemistry\nmatching and eligibility filters')
box(.735,.60,.245,.23,'Reference recovery','One bound geometry\nper CASF entry')
arrow((.27,.835),(.365,.745));arrow((.27,.555),(.365,.675));arrow((.62,.715),(.735,.715))
ax.text(.49,.485,'ChEMBL-count and 1,000-candidate targets\nPoseBusters-filtered CASF ensembles',ha='center',fontsize=8,color='#626C78')
box(.015,.13,.255,.20,'PLINDER','2024-06/v2\n51,280 ligand identities',color='#00866B')
box(.365,.13,.255,.20,'23-molecule panel','Physicochemical screen\n+ manual selection',color='#00866B')
box(.735,.13,.245,.20,'Coverage and precision','One ensemble compared\nwith all references',color='#00866B')
arrow((.27,.23),(.365,.23));arrow((.62,.23),(.735,.23))
fig.suptitle('Experimental references for protein-independent ensembles',x=.04,ha='left',y=.965,fontsize=12)
footer(fig,'Main panel: 94 CASF entries  •  Complementary panel: 23 molecules')
save(fig,'figure-1-datasets','**Figure S1. Experimental evaluation panels.** CASF-2016 is matched to ChEMBL3D by molecular identity and stereochemistry, followed by eligibility filtering, to obtain 94 core entries. The stored ChEMBL3D ensemble is the computed comparison resource. The separate 23-molecule PLINDER panel assesses several references per molecule. The larger 1,236-entry collection remains in Supporting Information. Arrows describe evaluation design, not a quantified exclusion funnel. For PLINDER 2024-06/v2, physicochemical screening reduced 51,280 identities to 25,392; ranking by distinct PDB entries and manual review yielded the 23-molecule panel (Methods; Appendix A). The supplied-pool comparison does not yet use the common CASF validity filter; reference-coordinate preprocessing remains to be documented.')

# 2. Count-matched geometric breadth and nearest-reference distance.
radii = pd.read_csv(TABLES/'clustering_radius_comparison.csv').set_index('method')
inputs.append(TABLES/'clustering_radius_comparison.csv')
fig,ax=plt.subplots(figsize=(7.2,4.8));fig.subplots_adjust(left=.12,right=.96,top=.86,bottom=.18)
axis(ax,'ChEMBL-count target · 94 CASF entries')
offsets={'chembl3d_gt_pb':(10,-20),'rdkit_random_raw':(10,3),'rdkit_random_minimized':(12,4),
 'torsion_raw':(-12,8),'torsion_minimized':(10,8),'loqi_raw':(-10,13),
 'torsional_diffusion_raw':(10,7),'mcf_drugs_l_raw':(12,9),'nextmol_dmt_l_raw':(12,-9),
 'flowr_raw':(-12,-4),Q:(-12,-17)}
plot_rows=[]
for key in KEYS:
 r=row(key,'chembl_count');x=float(radii.loc[r.method,'clusters_0p5']);y=float(r.best_rmsd)
 assert int(radii.loc[r.method,'n_defined']) == int(r.evaluated)
 plot_rows.append({'method':r.method,'label':STYLE[key][0],'clusters_0p5':x,'best_rmsd':y,'n_defined':int(r.evaluated)})
 point(ax,key,x,y)
 dx,dy=offsets[key]
 ax.annotate(STYLE[key][0],(x,y),xytext=(dx,dy),textcoords='offset points',fontsize=8,
  ha='right' if dx<0 else 'left',va='bottom',color=STYLE[key][1],
  arrowprops={'arrowstyle':'-','color':STYLE[key][1],'lw':.55,'alpha':.6})
ax.set(xlim=(15,62),ylim=(.45,.73),xlabel='Mean geometric clusters at 0.5 Å (more is better)',
 ylabel='Best RMSD (Å; lower is better)')
ax.set_xticks([20,30,40,50,60]);ax.set_yticks([.45,.50,.55,.60,.65,.70])
footer(fig,'Closest conformer per molecule · descriptive means · axes show the observed region')
pd.DataFrame(plot_rows).to_csv(OUT/'diversity-rmsd-points.csv',index=False)
save(fig,'figure-2-diversity-recovery','**Figure 1. Geometric diversity and proximity to the experimental bound conformation at the ChEMBL-count target.** Each point represents one evaluated pipeline or the stored ChEMBL3D-PB ensemble. The horizontal axis gives mean cluster count at 0.5 Å; the vertical axis gives Best RMSD, the per-molecule minimum over retained conformers averaged across entries with defined measurements. More clusters and lower Best RMSD place favorable ensembles toward the lower right. Means use 92 entries for NExT-Mol, 93 for MCF, and 94 for the others. Axes show the observed region for readability. Retained ensemble sizes differ despite matched candidate targets. These are descriptive means without uncertainty intervals; Best RMSD does not describe every generated conformer.')

# 3. Two endpoint comparisons, not interpolated sampling curves.
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,5.45),sharey=True)
fig.subplots_adjust(left=.235,right=.97,top=.86,bottom=.19,wspace=.25)
for i,key in enumerate(KEYS):
 r0=row(key,'chembl_count');r1=row(key,'fixed');color=STYLE[key][1]
 for ax,col,scale in [(a,'hit_0p75',100),(b,'clusters',1)]:
  x0=r0[col]*scale;x1=r1[col]*scale
  if key=='chembl3d_gt_pb':point(ax,key,x0,i,50)
  else:
   ax.plot([x0,x1],[i,i],color=color,lw=1.8,alpha=.6,zorder=2)
   ax.scatter(x0,i,s=33,facecolors='white',edgecolors=color,linewidth=1.3,zorder=3)
   ax.scatter(x1,i,s=37,color=color,edgecolors='white',linewidth=.5,zorder=4)
 for ax in (a,b):ax.axhspan(i-.42,i+.42,color='#F5F7F9' if i%2==0 else 'white',zorder=-1)
a.set_yticks(range(len(KEYS)),[STYLE[k][0] for k in KEYS]);a.set_ylim(len(KEYS)-.4,-.6)
a.set(xlim=(58,100),xlabel='Recovery at 0.75 Å (%)');a.set_xticks([60,70,80,90,100])
b.set(xlim=(0,125),xlabel='Mean clusters at 1.0 Å');b.set_xticks([0,40,80,120])
axis(a,'A   Experimental recovery','x');axis(b,'B   Geometric diversity','x')
for ax in(a,b):ax.tick_params(axis='y',length=0)
fig.legend(handles=[Line2D([],[],marker='o',color='#697581',mfc='white',ls='none',label='ChEMBL-count target'),
 Line2D([],[],marker='o',color='#697581',ls='none',label='1,000-candidate target'),
 Line2D([],[],marker='h',color='#333A43',ls='none',label='Stored ChEMBL3D-PB')],
 loc='lower center',bbox_to_anchor=(.53,.07),frameon=False,ncol=3,fontsize=7.5,handletextpad=.4,columnspacing=1.2)
footer(fig,'Lines join two evaluated targets; retained counts and computing costs are not matched.')
save(fig,'figure-3-sampling-budget','**Figure 2. Recovery and diversity at the two candidate targets.** Open circles denote the ChEMBL-count target and filled circles the 1,000-candidate target. Each row follows the same method between the two targets. The stored ChEMBL3D-PB ensemble is shown once as a hexagon. Recovery uses all 94 CASF entries; cluster means use defined measurements. Lines connect observed endpoints and do not represent random-subsampling curves or intermediate measurements. Targets precede PoseBusters filtering; retained counts are given in Table 3.')

# 4. Reconstruct strata from the original size/flexibility exports.
sizepath=ROOT/'docs/publication_tables_2026_09_28/size_strata.csv'
flexpath=ROOT/'docs/publication_tables_2026_09_24/flexibility.csv'
inputs.extend([sizepath,flexpath]);sizes=pd.read_csv(sizepath);flex=pd.read_csv(flexpath)
strata=[]
for key in KEYS:
 method=key if key=='chembl3d_gt_pb' else key+'_fixed'
 for desc,source,col,groups in [('size',sizes,'heavy_atoms',['<20','20-29','30-39','40+']),('rotors',flex,'rotors',['0-3','4-6','7-8','9+'])]:
  rows=source[source['set'].eq('core') & source.method.eq(method)]
  assert len(rows)==4,(method,desc)
  for r in rows.to_dict('records'):
   hits=r.get('hits',r['hit_075']*r['n']);assert np.isclose(hits,round(hits))
   strata.append(dict(key=key,descriptor=desc,group=r[col],n=int(r['n']),hits=int(round(hits)),hit=100*r['hit_075']))
  assert rows.n.sum()==94
strata=pd.DataFrame(strata);strata.to_csv(OUT/'strata.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(7.2,5.25),sharey=True)
fig.subplots_adjust(left=.235,right=.975,top=.85,bottom=.27,wspace=.18)
for j,(desc,groups,labels) in enumerate([('size',['<20','20-29','30-39','40+'],['<20','20–29','30–39','≥40']),('rotors',['0-3','4-6','7-8','9+'],['0–3','4–6','7–8','≥9'])]):
 ax=axes[j];d=strata[strata.descriptor.eq(desc)].set_index(['key','group'])
 matrix=np.array([[d.loc[(key,g),'hit'] for g in groups] for key in KEYS])
 mesh=ax.imshow(matrix,cmap='Blues',vmin=0,vmax=100,aspect='auto')
 for y,key in enumerate(KEYS):
  for x,g in enumerate(groups):
   r=d.loc[(key,g)];ax.text(x,y,f'{int(r.hits)}/{int(r.n)}',ha='center',va='center',fontsize=8,color='white' if r.hit>=65 else '#29323C')
 counts=[int(d.loc[(KEYS[0],g),'n']) for g in groups]
 ax.set_xticks(range(4),[f'{label}{"*" if n<5 else ""}\n(n = {n})' for label,n in zip(labels,counts)])
 ax.set_yticks(range(len(KEYS)),[STYLE[k][0] for k in KEYS]);ax.tick_params(which='both',length=0)
 ax.set_xticks(np.arange(-.5,4),minor=True);ax.set_yticks(np.arange(-.5,len(KEYS)),minor=True);ax.grid(which='minor',color='white',lw=1.2)
 ax.set_title('A   Heavy atoms' if j==0 else 'B   Rotatable bonds',loc='left',pad=13)
 for x,n in enumerate(counts):
  if n<5:ax.add_patch(Rectangle((x-.5,-.5),1,len(KEYS),fill=False,edgecolor='#CE6428',linestyle='--',linewidth=1.1,zorder=6,clip_on=False))
 for spine in ax.spines.values():spine.set_visible(False)
cax=fig.add_axes([.39,.14,.44,.024]);fig.colorbar(mesh,cax=cax,orientation='horizontal',label='Recovery at 0.75 Å (%)')
footer(fig,'Cells show recovered / total entries.  * Fewer than five entries: interpret descriptively.')
save(fig,'figure-4-size-flexibility','**Figure 3. Recovery by molecular size and flexibility.** Generated methods use the 1,000-candidate target; ChEMBL3D-PB uses its stored ensemble. Every cell gives recovered/total CASF entries, with the color indicating the corresponding percentage. Dashed outlines mark groups with fewer than five entries. Heavy-atom and rotatable-bond groups are analyzed separately and do not isolate independent effects of size and flexibility. Sparse groups cannot establish stable method rankings.')

# 5. Mean coverage/precision plus the paired 23-molecule comparison.
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,4.25));fig.subplots_adjust(left=.095,right=.975,top=.84,bottom=.23,wspace=.36)
axis(a,'A   Coverage and precision');axis(b,'B   Paired reference coverage')
d_offsets={'loqi_raw':(-5,11),Q:(-6,11),'torsional_diffusion_raw':(4,-17),'mcf_drugs_l_raw':(10,7),'nextmol_dmt_l_raw':(-10,23),'flowr_raw':(-10,-22)}
for key,label in DKEY.items():
 r=means.loc[label];x=r.cov_p_075*100;y=r.cov_r_075*100;point(a,key,x,y,60)
 dx,dy=d_offsets[key];a.annotate(STYLE[key][0],(x,y),xytext=(dx,dy),textcoords='offset points',fontsize=7.5,
  ha='right' if dx<0 else 'left',color=STYLE[key][1],arrowprops={'arrowstyle':'-','lw':.55,'color':STYLE[key][1]})
a.set(xlim=(15,76),ylim=(70,100),xlabel='Sample precision, COV-P (%)',ylabel='Reference coverage, COV-R (%)');a.set_xticks([20,40,60]);a.set_yticks([70,80,90,100])
q=drug[drug.label.eq(Q)].set_index('smiles');l=drug[drug.label.eq(DKEY['loqi_raw'])].set_index('smiles').reindex(q.index)
assert len(q)==len(l)==23 and l.cov_r_075.notna().all()
paired=pd.DataFrame({'name':q['name'],'LoQI':l.cov_r_075*100,'Qwen':q.cov_r_075*100});paired.to_csv(OUT/'paired_reference_coverage.csv',index=False)
points=paired.groupby(['LoQI','Qwen']).size().reset_index(name='n')
b.plot([0,100],[0,100],ls='--',color='#A1A9B1',lw=.8)
b.scatter(points.LoQI,points.Qwen,s=25+15*(points.n-1),color='#8A96A3',edgecolor='white',lw=.7,zorder=3)
for r in points[points.n>1].itertuples():
 b.annotate(f'{r.n} molecules',(r.LoQI,r.Qwen),xytext=(-12,-22),textcoords='offset points',ha='right',fontsize=7.5,
 arrowprops={'arrowstyle':'-','color':'#8A96A3','lw':.6})
# Names retained from the existing manuscript examples; not selected anew for advantage.
for name,offset in [('Imatinib',(8,7)),('Actinonin',(-8,12))]:
 r=paired[paired.name.str.lower().eq(name.lower())];assert len(r)==1;r=r.iloc[0]
 point(b,Q if r.Qwen>r.LoQI else 'loqi_raw',r.LoQI,r.Qwen,48)
 b.annotate(name,(r.LoQI,r.Qwen),xytext=offset,textcoords='offset points',fontsize=7.5,ha='right' if offset[0]<0 else 'left',
 arrowprops={'arrowstyle':'-','color':'#8A96A3','lw':.6})
b.set(xlim=(-5,106),ylim=(-5,106),xlabel='LoQI coverage (%)',ylabel='Qwen coverage (%)');b.set_xticks([0,50,100]);b.set_yticks([0,50,100])
footer(fig,'Exploratory · 23 molecules · supplied pools without a common PoseBusters filter')
save(fig,'figure-5-multiple-references','**Figure 4. Coverage of multiple bound references and precision relative to those observations.** (A) Molecule-averaged COV-R and COV-P for six learned generators, using the existing strict RMSD < 0.75 Å criterion. Axes are restricted to the observed range for readability. (B) Paired reference coverage for Qwen and LoQI across all 23 molecules; the diagonal denotes equal coverage. Coincident points are grouped and labeled, including 11 molecules with complete coverage by both methods. Imatinib and actinonin are the examples retained from the earlier draft. These are descriptive supplied-pool results with differing RMSD-failure conventions, not a comparison after common validity filtering. No uncertainty intervals are shown; paired uncertainty is documented in the archived findings.')

status='These figures use the archived corrected findings. Qwen is the initial 1.7B FSQ step-47,023 model, not the pending benchmark-excluded model. The corrected energy figure and dataset-release figures remain pending.'
gallery='# Results figures\n\n'+status+'\n\n'
for name,caption in figures:
 gallery+=f'![{name}]({OUT/name}.png)\n\n{caption}\n\n[PDF]({OUT/name}.pdf) · [Editable SVG]({OUT/name}.svg)\n\n'
(OUT/'gallery.md').write_text(gallery)
(OUT/'provenance.json').write_text(json.dumps({'prepared_on':'2026-10-05','qwen_checkpoint':Q,'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'style':STYLE,'figure_files':[n for n,_ in figures],'notes':status},indent=2)+'\n')
# Place figures at the existing placeholders. Keep inserted figures outside table markers.
p=ROOT/'docs/results_working.md';s=p.read_text();discussion=s.split('# Discussion\n',1)[1]
methods_path=ROOT/'docs/appendix_working.md'
methods_text=methods_path.read_text() if methods_path.exists() else None
prefixes=['[Here goes the dataset construction diagram,','[Here goes the cluster-count versus recovery figure',
 '[Here goes the figure connecting recovery at the two candidate targets.',
 '[Here goes the figure showing recovery at the 1,000-candidate target by',
 '[Here goes the coverage-versus-precision figure with representative']
for i,((name,caption),prefix) in enumerate(zip(figures,prefixes),1):
 block=f'<!-- results-figure-{i}:start -->\n![Figure {"S1" if i == 1 else i-1}]({OUT/name}.png)\n\n{caption}\n\n[PDF]({OUT/name}.pdf) · [Editable SVG]({OUT/name}.svg)\n<!-- results-figure-{i}:end -->'
 pattern=rf'<!-- results-figure-{i}:start -->.*?<!-- results-figure-{i}:end -->'
 if i == 1 and methods_text is not None:
  methods_text,n=re.subn(pattern,lambda _:block,methods_text,flags=re.S)
  assert n==1,'Dataset figure marker missing from appendix draft'
  continue
 if re.search(pattern,s,re.S):s=re.sub(pattern,lambda _:block,s,flags=re.S)
 elif i in [3,5]:
  s=re.sub(re.escape(prefix)+r'[^\n]*\]', '', s)
  marker=f'<!-- results-table-{4 if i==3 else 5}:end -->'
  assert marker in s;s=s.replace(marker,marker+'\n\n'+block,1)
 else:
  s,n=re.subn(re.escape(prefix)+r'[^\n]*\]',lambda _:block,s,count=1);assert n==1,i
for prefix in (prefixes[2],prefixes[4]):
 s=re.sub(re.escape(prefix)+r'[^\n]*\]', '', s)
assert s.split('# Discussion\n',1)[1]==discussion
p.write_text(s)
if methods_text is not None:methods_path.write_text(methods_text)
print(f'Created {len(figures)} figures in PDF, SVG, and 350-dpi PNG; gallery and draft updated.')
