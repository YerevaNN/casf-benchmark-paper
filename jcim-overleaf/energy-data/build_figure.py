"""Render the energy comparison from the recorded common-hydrogen rescoring."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'docs/results_figures'
rows=list(csv.DictReader((Path(__file__).parent/'summary.csv').open()))
style=json.loads((OUT/'provenance.json').read_text())['style']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.titleweight':'medium','axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.edgecolor':'#8B929A','axes.labelcolor':'#29323C','text.color':'#29323C',
 'xtick.color':'#535D68','ytick.color':'#535D68','axes.linewidth':.65,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
 'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,5.4),sharey=True)
fig.subplots_adjust(left=.235,right=.965,top=.85,bottom=.17,wspace=.25)
y=np.arange(len(rows))
for ax,key,title,xlabel in [(a,'mean_clusters_0p5','A   Geometric diversity','Mean clusters at 0.5 Å'),
 (b,'median_energy_std','B   Energy dispersion','Median energy SD (kcal/mol)')]:
 vals=[float(r[key]) for r in rows]
 colors=[style[r['method'].removesuffix('_chembl_count')][1] for r in rows]
 ax.barh(y,vals,color=colors,height=.58,zorder=3)
 ax.set_xlim(0,max(vals)*1.2)
 ax.set_title(title,loc='left',pad=14)
 ax.set_xlabel(xlabel,labelpad=10)
 ax.grid(axis='x',color='#E7EBEE',linewidth=.6,zorder=0)
 ax.set_axisbelow(True)
 ax.tick_params(axis='y',length=0)
 ax.axhline(len(rows)-1.5,color='#7B838D',lw=.7,ls=(0,(4,3)))
 for i,v in enumerate(vals):ax.text(v+max(vals)*.025,i,f'{v:.1f}',va='center',fontsize=8)
a.set_yticks(y,labels=[r['label'] for r in rows]);a.invert_yaxis()
b.tick_params(labelleft=False)
fig.suptitle('ChEMBL-count target · geometry and energy spread',x=.235,ha='left',y=.96,fontsize=11)
fig.text(.6,.052,'Common MMFF94s protocol · unoptimized hydrogens · original heavy-atom coordinates',ha='center',fontsize=7.3,color='#626C78')
for ext in ['pdf','svg','png']:fig.savefig(OUT/('figure-energy-dispersion.'+ext),dpi=350)
plt.close(fig)
svg=OUT/'figure-energy-dispersion.svg';svg.write_text('\n'.join(x.rstrip() for x in svg.read_text().splitlines())+'\n')
CAPTION='**Figure 2. Geometric diversity and energy dispersion at the ChEMBL-count target.** (A) Mean geometric cluster count at 0.5 Å. (B) Median across molecules of the population standard deviation of conformer energies within each retained ensemble, in kcal/mol. Energies were recalculated with MMFF94s after rebuilding hydrogen coordinates without optimizing any atom. Each method uses the same measurable entries in both panels: 92 for NExT-Mol, 93 for MCF, and 94 for the others. The dashed line separates the stored ChEMBL3D-PB ensemble. Smaller energy SD indicates a narrower distribution; it does not establish lower absolute energy or the absence of high-energy conformers. Means, paired energy differences, and calculation completeness are reported in Supporting Table S9.'
(OUT/'figure-energy-dispersion-caption.md').write_text(CAPTION+'\n')
# Preserve the manuscript's placement when either renderer is rerun.
import re
name='figure-energy-dispersion'
block=f'<!-- results-figure-energy:start -->\n![Figure 2]({OUT/name}.png)\n\n{CAPTION}\n\n[PDF]({OUT/name}.pdf) · [Editable SVG]({OUT/name}.svg)\n<!-- results-figure-energy:end -->'
draft=ROOT/'docs/results_working.md'
s=draft.read_text()
pattern=r'<!-- results-figure-energy:start -->.*?<!-- results-figure-energy:end -->'
if re.search(pattern,s,re.S):s=re.sub(pattern,lambda _:block,s,flags=re.S)
else:
    anchor='Constructing a large conformational resource also requires'
    assert anchor in s
    s=s.replace(anchor,block+'\n\n'+anchor,1)
draft.write_text(s)
gallery=OUT/'gallery.md';s=gallery.read_text()
if re.search(pattern,s,re.S):s=re.sub(pattern,lambda _:block,s,flags=re.S)
else:
    anchor='![figure-3-sampling-budget]'
    assert anchor in s
    s=s.replace(anchor,block+'\n\n'+anchor,1)
gallery.write_text(s)
print('Rendered energy figure and updated working-draft displays.')
