"""Validate selected SI against current main-text data and LaTeX inclusion order."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
PAPER=HERE.parent
manifest=json.loads((HERE/'provenance.json').read_text())
for section in ['inputs','outputs']:
    for name,digest in manifest[section].items():
        assert hashlib.sha256((PAPER/name).read_bytes()).hexdigest()==digest,name

changes=pd.read_csv(HERE/'recovery-gains-losses.csv')
assert len(changes)==20
np.testing.assert_array_equal(changes[['gained','lost','shared_hits','shared_misses']].sum(axis=1),changes.n)
assert changes.loc[(changes.cohort=='ref')&(changes.name=='Qwen 1.7B FSQ'),['gained','lost']].iloc[0].tolist()==[155,15]

counts=pd.read_csv(PAPER/'results-data/three-dataset-tables/table3.csv')
styles=json.loads((PAPER/'sampling-data/palette.json').read_text())
for cohort,path in [('core','figure-data/candidates/core-entries.csv'),('ref','results-data/ref-panel/per-entry.csv')]:
    source=pd.read_csv(PAPER/path)
    for key,(name,*_) in styles.items():
        for target in (['chembl_count'] if key=='chembl3d_gt_pb' else ['chembl_count','fixed']):
            method=key if key=='chembl3d_gt_pb' else key+'_'+target
            rows=source[source.method==method]
            assert len(rows)==(94 if cohort=='core' else 1044)
            expected=counts[(counts.cohort==cohort)&(counts.method==name)&(counts.budget==target)].iloc[0]
            np.testing.assert_allclose(rows.post_pb_confs.fillna(0).mean(),expected.mean_retained,atol=1e-10)
            tex=(PAPER/'tables/supporting-ensemble-counts.tex').read_text()
            line=next(line for line in tex.splitlines() if line.startswith(name+' &'))
            index=(0 if cohort=='core' else 2)+(0 if target=='chembl_count' else 1)
            assert line.split(' & ')[index+1].replace(r' \\','')==f'{expected.mean_retained:.1f}'

pl=pd.read_csv(HERE/'plinder-by-molecule.csv')
assert len(pl)==242 and pl.label.nunique()==11 and pl.smiles.nunique()==22
assert not pl.label.isin(['torsion_raw_druglike','torsion_minimized_22_druglike']).any()
q=pd.read_csv(HERE/'qwen-validity-source.csv')
assert len(q)==88 and q.model.nunique()==4
matched=pl[pl.label=='qwen_1p7b_v2_coordtoken_step41054'].merge(q[q.model=='fsq17'],on='smiles',validate='one_to_one',suffixes=('_main','_comparison'))
for column in ['n_generated','n_pb_pass']:
    np.testing.assert_array_equal(matched[column+'_main'],matched[column+'_comparison'])

# Expand includes in place to check the automatic numbering in the separate SI.
def expand(path):
    text=re.sub(r'(?<!\\)%[^\n]*','',path.read_text())
    def include(match):
        name=match[1].replace(r'\manuscriptroot','').strip()
        return expand(PAPER/(name if name.endswith('.tex') else name+'.tex'))
    return re.sub(r'\\input\{([^}]+)\}',include,text)
si=expand(PAPER/'supporting-information.tex')
expected_tables=['supporting-ensemble-counts','plinder23-screen','plinder23-statistics','plinder23-molecules','generators','energy-summary']
tables=re.findall(r'\\label\{tab:([^}]+)\}',si)
# Exact IDs are retained from the existing source files.
assert tables==expected_tables,tables
figures=re.findall(r'\\label\{fig:([^}]+)\}',si)
assert len(figures)==6,figures
assert figures[1:]==['qwen-molecule-effects','supporting-gains-losses','supporting-composition','supporting-plinder-coverage','supporting-plinder-validity'],figures
assert 'figures/candidates/' not in si
for forbidden in ['Mean median','tab:supporting-druglike','tab:size','tab:supporting-rmsd']:
    assert forbidden not in si
main=expand(PAPER/'acs-template.tex')
assert 'figures/candidates/' not in main and 'Additional figure candidates for review' not in main
for source in [main,si]:
    for number in re.findall(r'(?:Supporting )?(?:Table|Figure)s?\s*~?\s*S(\d+)',source):
        assert 1<=int(number)<=6,number
    for name in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',source):
        assert (PAPER/name).is_file(),name
print('Verified current source/artwork hashes, 20 paired recovery comparisons, retained counts, 242 PLINDER and 88 v2 validity records.')
print('Current SI: six tables and six figures; no historical review gallery or removed tables included.')
print('SI table order:',tables)
print('SI figure order:',figures)
