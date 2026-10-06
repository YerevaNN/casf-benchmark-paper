"""Render the approved threshold curves from archived values."""
from pathlib import Path
import json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
STYLE=json.loads((OUT/'palette.json').read_text())
curves=pd.read_csv(OUT/'threshold_curves.csv')
METHODS={
'qwen_1p7b_fsq_bigdata_step47023':'qwen_1p7b_fsq_bigdata_step47023',
'loqi_druglike':'loqi_raw',
'torsional_diffusion_druglike':'torsional_diffusion_raw',
'mcf_drugs_l_druglike':'mcf_drugs_l_raw',
'flowr_druglike':'flowr_raw',
'nextmol_dmt_l_druglike':'nextmol_dmt_l_raw',
}
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
fig.text(.5,.02,'Strict RMSD < threshold. Supplied ensembles without common validity filtering.',ha='center',fontsize=8,color='#626C78')
for ext in ['png','pdf','svg']:fig.savefig(OUT/('figure-5-threshold-coverage.'+ext),dpi=250)
plt.close(fig)
