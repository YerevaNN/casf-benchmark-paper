"""Check current paper displays against the archived PLINDER-22 measurements."""
from pathlib import Path
import csv
import io
import json
import re
import subprocess
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
REPO = PAPER.parent
methods = json.loads((HERE / 'methods.json').read_text())
frame = pd.read_csv(HERE / 'per-molecule.csv')
summary = pd.read_csv(HERE / 'summary.csv').set_index('label')
means = pd.read_csv(PAPER / 'results-data/druglike_means.csv').set_index('label')
curves = pd.read_csv(HERE / 'threshold-curves.csv')
inventory = pd.read_csv(HERE / 'reference-inventory.csv')
records = pd.read_csv(HERE / 'reference-records.csv')
assert len(frame) == 22 * len(methods) and len(inventory) == 22 and len(records) == 1067
assert inventory.num_true_confs.sum() == 1067
assert records.pdb_id.nunique() == 552 and records.system_id.nunique() == 1045
assert not records.duplicated(['ccd','conf_index']).any()
assert 'BO2' not in set(records.ccd)
assert set(frame.smiles) == set(inventory.smiles)
entries = json.loads((HERE / 'entries.json').read_text())
with np.load(HERE / 'nearest-distances.npz', allow_pickle=False) as archive:
    for method, label in methods.items():
        group = frame.loc[frame.label.eq(method)].set_index('smiles')
        assert len(group) == 22 and group.num_true_confs.sum() == 1067
        assert group.conformer_count.eq(group.n_generated).all()
        plotted = curves.loc[curves.method.eq(method)].sort_values('threshold')
        thresholds = plotted.threshold.to_numpy()
        recall, precision = [], []
        for e in [e for e in entries if e['method'] == method]:
            row = group.loc[e['smiles']]
            r, p = archive[e['key']+'_reference'], archive[e['key']+'_generated']
            assert np.isfinite(r).all() and np.isfinite(p).all()
            assert len(r) == row.num_true_confs and len(p) == row.num_gen_confs
            np.testing.assert_allclose([r.mean(),p.mean()], [row.mat_r,row.mat_p], rtol=0,atol=1e-12)
            recall.append((r[:,None] <= thresholds).mean(axis=0))
            precision.append((p[:,None] <= thresholds).mean(axis=0))
        np.testing.assert_allclose(np.mean(recall,axis=0), plotted.recall, rtol=0,atol=1e-12)
        np.testing.assert_allclose(np.mean(precision,axis=0), plotted.precision, rtol=0,atol=1e-12)
        for column in ['recall','precision']:
            assert (np.diff(plotted[column]) >= -1e-12).all()
        for col,source in [('num_gen_confs','num_gen_confs'),('cov_r_075','cov_r_0p75'),
                           ('cov_p_075','cov_p_0p75'),('mat_r','mat_r'),('mat_p','mat_p')]:
            np.testing.assert_allclose(means.loc[method,col],group[source].mean(),rtol=0,atol=1e-12)
        np.testing.assert_allclose(means.loc[method,'pooled_pb_pass'],group.n_pb_pass.sum()/group.n_generated.sum(),rtol=0,atol=1e-12)
        for number, columns in [(1,{'clusters_0p5':'greedy_clusters_0p5','clusters_1p0':'greedy_clusters_1p0','largest_percent':'largest_cluster_fraction_1p0'}),
                                (2,{'recovery_0p5':'cov_r_0p5','recovery_0p75':'cov_r_0p75','distance':'mat_r'})]:
            table = pd.read_csv(PAPER / f'results-data/three-dataset-tables/table{number}.csv')
            row = table.loc[table.cohort.eq('plinder') & table.method.eq(label)].iloc[0]
            assert row.n_total == 22
            for dest,source in columns.items():
                scale = 100 if dest == 'largest_percent' or dest.startswith('recovery') else 1
                np.testing.assert_allclose(row[dest],scale*group[source].mean(),rtol=0,atol=1e-12)

# Core values remain the approved snapshot. Ref is checked against its selected cohort.
for number in [1,2,3]:
    path = f'jcim-overleaf/results-data/three-dataset-tables/table{number}.csv'
    prior = subprocess.check_output(['git','show',f'c797287:{path}'],cwd=REPO,text=True)
    old = list(csv.DictReader(io.StringIO(prior)))
    new = list(csv.DictReader((REPO/path).open()))
    assert [r for r in old if r['cohort'] == 'core'] == [r for r in new if r['cohort'] == 'core']
# Verify rendered table cell text against CSVs (including all missing cells).
for number,name,columns in [(1,'diversity',['clusters_0p5','clusters_1p0','largest_percent']),
                            (2,'count-recovery',['recovery_0p5','recovery_0p75','distance'])]:
    data = pd.read_csv(PAPER / f'results-data/three-dataset-tables/table{number}.csv')
    tex = (PAPER / f'tables/{name}.tex').read_text()
    for method in data.method.unique():
        line = next(line for line in tex.splitlines() if line.startswith(method+' &') or line.startswith(method+'$'))
        cells = line.split(' & ')[1:]
        for ci,cohort in enumerate(['plinder','core','ref']):
            row = data.loc[data.method.eq(method) & data.cohort.eq(cohort)].iloc[0]
            for j,column in enumerate(columns):
                expected = '---' if pd.isna(row[column]) else f'{row[column]:.{3 if column=="distance" else 1}f}'
                cell = cells[ci*3+j].replace(r' \\', '').strip()
                cell = re.sub(r'\\(?:textbf|underline)\{([^}]+)\}', r'\1', cell)
                assert cell == expected, (name,method,cohort,column,cell,expected)
        if method == 'ChEMBL3D-PB':
            assert r'\textbf' not in line and r'\underline' not in line
for script,names in [('results-data/three-dataset-tables/build_tables.py',['diversity','count-recovery','recovery']),
                     ('results-data/build_druglike_table.py',['druglike','supporting-druglike']),
                     ('plinder22-data/build_inventory.py',['plinder23-statistics','plinder23-molecules'])]:
    before = {n:(PAPER/f'tables/{n}.tex').read_bytes() for n in names}
    import sys
    subprocess.run([sys.executable,str(PAPER/script)],check=True,capture_output=True)
    assert all((PAPER/f'tables/{n}.tex').read_bytes()==value for n,value in before.items())

# Source checks for the main paper and separate supporting document; no compilation.
bibkeys=set(re.findall(r'@\w+\s*\{\s*([^,]+)',(PAPER/'acs-template.bib').read_text()))
for root in ['acs-template.tex','supporting-information.tex']:
    seen=set()
    def visit(path):
        if path in seen: return ''
        seen.add(path)
        text=re.sub(r'(?<!\\)%[^\n]*','',path.read_text())
        # Brace balance, ignoring escaped braces.
        depth=0
        for c in re.findall(r'(?<!\\)[{}]',text):
            depth += 1 if c=='{' else -1
            assert depth>=0,(path,'closing brace')
        assert depth==0,(path,'brace balance',depth)
        stack=[]
        for action,env in re.findall(r'\\(begin|end)\{([^}]+)\}',text):
            if action=='begin':stack.append(env)
            else: assert stack.pop()==env,(path,env)
        assert not stack,(path,stack)
        extra=[]
        for target in re.findall(r'\\input\{([^}]+)\}',text):
            target=target.replace(r'\manuscriptroot','').strip()
            if not target.endswith('.tex'):target+='.tex'
            extra.append(visit(PAPER/target))
        for target in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',text):
            assert (PAPER/target).is_file(),(path,target)
        return text+'\n'+'\n'.join(extra)
    text=visit(PAPER/root)
    labels=re.findall(r'\\label\{([^}]+)\}',text)
    assert len(labels)==len(set(labels)),(root,'duplicate labels')
    refs=set(re.findall(r'\\(?:ref|eqref|autoref)\{([^}]+)\}',text))
    assert refs.issubset(labels),(root,refs-set(labels))
    citations={key.strip() for group in re.findall(r'\\cite\w*\{([^}]+)\}',text) for key in group.split(',')}
    assert citations.issubset(bibkeys),(root,citations-bibkeys)
    print(f'{root}: {len(seen)} source files, {len(labels)} labels, {len(refs)} references, {len(citations)} citations checked.')
print('Verified 220 PLINDER method–molecule rows, all curves, table cells, inventory and unchanged CASF core data.')
