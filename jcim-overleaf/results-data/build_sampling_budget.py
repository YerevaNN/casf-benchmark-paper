"""Export the two-budget comparison from the pinned core94 measurements.
Run from the analysis checkout; no conformers are regenerated or reclustered.
"""
from pathlib import Path
import hashlib,json,sqlite3
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
DB=ROOT/'data/results/casf_analysis_dashboard.sqlite'
EXPECTED='bfd6093b269333f00864e7040bbac47b82dbd4c0f2b682659af94774cb7ab116'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(DB)==EXPECTED
selected=pd.read_csv(OUT/'casf_selected.csv')
with sqlite3.connect(f'file:{DB}?mode=ro',uri=True) as con:
    per=pd.read_sql_query('SELECT method,mol_id,casf_best_rmsd,greedy_clusters_0p5,post_pb_confs FROM per_ligand_long WHERE ligand_set="core"',con)
per=per[per.method.isin(selected.method)]
assert not per.duplicated(['method','mol_id']).any()
records=[]
for r in selected.itertuples():
    d=per[per.method.eq(r.method)]
    assert len(d)==94
    assert d.casf_best_rmsd.notna().equals(d.greedy_clusters_0p5.notna())
    best=d.casf_best_rmsd.mean();hit=d.casf_best_rmsd.le(.75).mean()
    assert np.isclose(best,r.best_rmsd) and np.isclose(hit,r.hit_0p75)
    assert np.isclose(d.post_pb_confs.mean(),r.mean_valid)
    records.append(dict(method=r.method,label=r.display_label,tier=r.tier,
        best_rmsd=best,hit_0p75=hit,clusters_0p5=d.greedy_clusters_0p5.mean(),
        n_defined=int(d.casf_best_rmsd.notna().sum()),mean_retained=d.post_pb_confs.mean()))
s=pd.DataFrame(records)
s.to_csv(OUT/'sampling_budget_comparison.csv',index=False)
per.to_csv(OUT/'sampling_budget_per_entry.csv',index=False)
# Ordering matches Table 2: lowest Best RMSD immediately above the stored reference.
fixed=s[s.tier.eq('fixed')].sort_values('best_rmsd',ascending=False)
ordered=[]
for r in fixed.itertuples():
    count=s[s.method.eq(r.method.removesuffix('_fixed')+'_chembl_count')].iloc[0]
    ordered.append((r.label,r,count))
ref=s[s.method.eq('chembl3d_gt_pb')].iloc[0]
ordered.append((ref.label,None,ref))
metrics=['best_rmsd','hit_0p75','clusters_0p5']
# Rank unrounded values, including the stored reference only where it is reported.
ranks={}
for tier in ['fixed','chembl_count']:
    rows=s[s.tier.eq(tier) | (s.tier.eq('reference') if tier=='chembl_count' else False)]
    for key in metrics:ranks[tier,key]=sorted(set(rows[key]),reverse=key!='best_rmsd')[:2]
def cell(r,key,tier,tex=False):
    if r is None:return '---' if tex else '—'
    value=float(r[key] if isinstance(r,pd.Series) else getattr(r,key))
    value_text=f'{value:.3f}' if key=='best_rmsd' else f'{value*(100 if key=="hit_0p75" else 1):.1f}'
    best=ranks[tier,key]
    if value==best[0]:return r'\textbf{'+value_text+'}' if tex else '**'+value_text+'**'
    if value==best[1]:return r'\underline{'+value_text+'}' if tex else '<u>'+value_text+'</u>'
    return value_text
note='Hit@0.75 uses all 94 entries, with missing outputs counted as failures. Best RMSD is the mean of the per-entry minimum RMSD; cluster counts are also averaged over entries with defined measurements. At the 1,000-candidate target, these means use 93 entries for MCF and NExT-Mol and 94 for the others; at the ChEMBL-count target, they use 93 for MCF, 92 for NExT-Mol, and 94 for the others. Generators are sorted by decreasing unrounded Best RMSD at 1,000 candidates. Bold and underlining identify the best and second-best distinct values in each column, including ties and the stored reference where reported. ChEMBL3D-PB is the same stored ensemble and has no 1,000-candidate result. Candidate targets precede filtering and do not match retained counts or computational cost.'
md=['**Table 3. Recovery and geometric diversity at the two candidate targets.**','',
'| Method | 1,000: Best RMSD (Å) ↓ | 1,000: Hit@0.75 (%) ↑ | 1,000: clusters at 0.5 Å ↑ | ChEMBL-count: Best RMSD (Å) ↓ | ChEMBL-count: Hit@0.75 (%) ↑ | ChEMBL-count: clusters at 0.5 Å ↑ |',
'| --- | --- | --- | --- | --- | --- | --- |']
tex=[r'% Source: results-data/sampling_budget_comparison.csv.',r'\begin{table}[htbp]',r'\centering\footnotesize',r'\setlength{\tabcolsep}{4pt}',r'\renewcommand{\arraystretch}{1.15}',r'\caption{Recovery and geometric diversity at the two candidate targets on the 94-entry core panel.}',r'\label{tab:recovery}',r'\begin{tabular}{@{}lrrr@{\hspace{14pt}}rrr@{}}',r'\toprule',r'& \multicolumn{3}{c}{1,000-candidate target} & \multicolumn{3}{c}{ChEMBL-count target} \\',r'\cmidrule(lr){2-4}\cmidrule(l){5-7}',r'Method & \shortstack{Best RMSD\\(\AA{}) $\downarrow$} & \shortstack{Hit@0.75\\(\%) $\uparrow$} & \shortstack{Clusters\\at 0.5 \AA{} $\uparrow$} & \shortstack{Best RMSD\\(\AA{}) $\downarrow$} & \shortstack{Hit@0.75\\(\%) $\uparrow$} & \shortstack{Clusters\\at 0.5 \AA{} $\uparrow$} \\',r'\midrule']
for label,f,c in ordered:
    if f is None:
        md.append('| ' + ' | '.join(['---']*7)+' |');tex.append(r'\hdashline')
    md.append('| '+' | '.join([label]+[cell(r,k,t) for t,r in [('fixed',f),('chembl_count',c)] for k in metrics])+' |')
    tex.append(' & '.join([label]+[cell(r,k,t,True) for t,r in [('fixed',f),('chembl_count',c)] for k in metrics])+r' \\')
md.extend(['',note,'','Source: [budget comparison](results_tables/sampling_budget_comparison.csv); [per-entry records](results_tables/sampling_budget_per_entry.csv).'])
tex.extend([r'\bottomrule',r'\end{tabular}',r'\par\smallskip',r'\begin{minipage}{\linewidth}',note.replace('Å',r'\AA{}'),r'\end{minipage}',r'\end{table}'])
(OUT/'sampling_budget_table.tex').write_text('\n'.join(tex)+'\n')
p=ROOT/'docs/results_working.md';text=p.read_text();a=text.index('<!-- results-table-4:start -->');b=text.index('<!-- results-table-4:end -->',a)
p.write_text(text[:a]+'<!-- results-table-4:start -->\n'+'\n'.join(md)+'\n'+text[b:])
(OUT/'sampling-budget-provenance.json').write_text(json.dumps(dict(database=str(DB),database_sha256=EXPECTED,selected_sha256=sha(OUT/'casf_selected.csv'),script_sha256=sha(Path(__file__)),cluster_radius_angstrom=.5,recovery_cutoff_angstrom=.75,sort='Decreasing unrounded 1,000-candidate Best RMSD; stored reference last',sources={p.name:sha(p) for p in [OUT/'sampling_budget_comparison.csv',OUT/'sampling_budget_per_entry.csv']}),indent=2)+'\n')
print(s[['label','tier','best_rmsd','hit_0p75','clusters_0p5']].to_string(index=False))
