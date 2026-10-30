"""Evaluation-panel diagram in the manuscript's existing box-and-arrow style."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
PAPER = Path(__file__).resolve().parent.parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
fig, ax = plt.subplots(figsize=(10,6.3))
fig.subplots_adjust(left=.02,right=.98,top=.89,bottom=.08)
ax.set(xlim=(0,1),ylim=(0,1)); ax.axis('off')
def box(x,y,w,h,title,body,color):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012,rounding_size=0.015',
        linewidth=1.2,edgecolor=color,facecolor='#F5F8FA'))
    ax.text(x+w/2,y+h*.74,title,ha='center',va='center',fontsize=11,color=color)
    ax.text(x+w/2,y+h*.35,body,ha='center',va='center',fontsize=9,color='#29323C',linespacing=1.4)
def arrow(start,end):
    ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'-|>','lw':1.2,'color':'#8B929A'})
green,blue='#00866B','#0072B2'
box(.02,.67,.25,.27,'PLINDER','2024-06/v2 selection\n51,280 → 25,392 identities',green)
box(.375,.67,.25,.27,'22-molecule panel','Manual selection + cleaning\nCovalent and stereo exclusions',green)
box(.73,.67,.25,.27,'Coverage and precision','1,067 deduplicated references\nMultiple conformations / molecule',green)
arrow((.28,.805),(.365,.805));arrow((.635,.805),(.72,.805))
ax.text(.5,.59,'Full supplied ensembles · equal molecule weights · validity assessed separately',ha='center',color='#626C78',fontsize=8.5)
box(.02,.15,.25,.27,'CASF + ChEMBL3D','Identity and stereochemistry\nmatching; eligibility filters',blue)
box(.375,.15,.25,.27,'Core: 94 · Ref: 1,044','Stored ChEMBL3D ensembles\nprovide the comparison',blue)
box(.73,.15,.25,.27,'Reference recovery','One bound geometry\nper CASF entry',blue)
arrow((.28,.285),(.365,.285));arrow((.635,.285),(.72,.285))
ax.text(.5,.06,'ChEMBL-count and 1,000-candidate targets · PoseBusters-filtered ensembles',ha='center',color='#626C78',fontsize=8.5)
fig.suptitle('Experimental references for protein-independent ensembles',x=.04,ha='left',y=.98,fontsize=13,color='#29323C')
for extension in ['pdf','png','svg']:
    fig.savefig(PAPER/'figures'/f'experimental-panels.{extension}',dpi=250)
plt.close(fig)

# Keep generated vector files clean in source control.
svg = PAPER / 'figures/experimental-panels.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
