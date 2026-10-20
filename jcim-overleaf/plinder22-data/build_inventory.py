"""Render PLINDER SI inventory from cleaned counts and archived identity descriptors."""
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
original = pd.read_csv(PAPER / 'method-data/plinder23-molecules.csv')
clean = pd.read_csv(HERE / 'reference-inventory.csv')
old = original.set_index('ccd').loc[clean.ccd].reset_index()
assert set(old.ccd) == set(clean.ccd) and 'BO2' not in set(clean.ccd)
# Preserve selection frequencies separately from counts in the cleaned SDFs.
old = old.rename(columns={'n_pdbs':'selection_n_pdbs', 'n_systems':'selection_n_systems',
                          'num_true_confs':'original_num_true_confs'})
old = old.drop(columns='smiles').merge(clean.drop(columns='name'), on='ccd', validate='one_to_one')
old = old.sort_values('rank')
old['selection_rank'] = old['rank']
old['rank'] = range(1,23)
old.to_csv(PAPER / 'method-data/plinder22-molecules.csv', index=False)
f = PAPER / 'tables/plinder23-molecules.tex'
s = f.read_text(); start = s.index('1 & Staurosporine'); end = s.index('\\end{longtable}', start)
rows = []
for r in old.itertuples():
    rows.append(f'{r.rank} & {r.name} & {r.ccd} & {r.n_pdbs} & {r.n_systems} & {r.num_true_confs} & {r.mw:.1f} & {r.heavy_atoms} & {r.rotb}' + r'\\')
s = s[:start] + '\n'.join(rows) + '\n' + s[end:]
s = s.replace('The 23 selected molecules in saved PDB-frequency order. CCD denotes the Chemical Component Dictionary code; reference counts describe the evaluated coordinate collection.', 'The 22 retained molecules in their original selection order. CCD denotes the Chemical Component Dictionary code. PDB, system, and reference counts describe the cleaned coordinate collection; original selection frequencies are archived separately.')
f.write_text(s)
f = PAPER / 'tables/plinder23-statistics.tex'; s = f.read_text(); start = s.index('Distinct PDB entries per identity'); end = s.index('\\end{longtable}',start)
rows=[]
for label,column,digits in [('Distinct PDB entries per identity','n_pdbs',0),('PLINDER systems per identity','n_systems',0),('Evaluated reference conformers per identity','num_true_confs',0),('Molecular weight (Da)','mw',1),('Heavy atoms','heavy_atoms',0),('Rotatable bonds','rotb',0),('Crippen cLogP','clogp',2),(r'TPSA (\AA{}\textsuperscript{2})','tpsa',1)]:
    values = [old[column].min(),old[column].median(),old[column].max()]
    def fmt(x):
        return str(int(x)) if digits==0 and x==int(x) else f'{x:.{max(digits,1)}f}'
    rows.append(label+' & '+' & '.join(map(fmt,values))+r'\\')
s=s[:start]+'\n'.join(rows)+'\n'+s[end:]
f.write_text(s.replace('Plinder-23 panel characteristics','Cleaned PLINDER-22 panel characteristics'))
print('Rendered the 22-molecule inventory and panel statistics.')
