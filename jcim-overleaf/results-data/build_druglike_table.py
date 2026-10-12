"""Render Table 4 from the archived per-molecule records. No scoring is run."""
from pathlib import Path
import csv
import math
import re
from statistics import mean

DATA = Path(__file__).resolve().parent
PAPER = DATA.parent
METHODS = [
    ('rdkit_random_raw_druglike', 'RDKit raw'),
    ('torsion_raw_druglike', 'Torsion raw'),
    ('loqi_druglike', 'LoQI'),
    ('torsional_diffusion_druglike', 'Torsional Diffusion'),
    ('mcf_drugs_l_druglike', 'MCF drugs-L'),
    ('nextmol_dmt_l_druglike', 'NExT-Mol DMT-L'),
    ('flowr_druglike', 'FlowR'),
    ('qwen_1p7b_fsq_bigdata_step47023', 'Qwen 1.7B FSQ'),
]
COLUMNS = ['num_gen_confs', 'cov_r_075', 'cov_p_075', 'mat_r', 'mat_p', 'pooled_pb_pass']
LOWER_IS_BETTER = {'mat_r', 'mat_p'}
with (DATA / 'druglike_selected.csv').open() as handle:
    records = list(csv.DictReader(handle))
assert {r['label'] for r in records} == {m for m, _ in METHODS}
rows = []
identities = None
for method, label in METHODS:
    selected = [r for r in records if r['label'] == method]
    current = {r['smiles'] for r in selected}
    assert len(selected) == len(current) == 23
    if identities is None:
        identities = current
    assert current == identities
    assert sum(int(r['num_true_confs']) for r in selected) == 2450
    values = {k: mean(float(r[k]) for r in selected) for k in COLUMNS[:-1]}
    generated = sum(int(r['n_generated']) for r in selected)
    passed = sum(int(r['n_pb_pass']) for r in selected)
    assert generated > 0 and 0 <= passed <= generated
    values['pooled_pb_pass'] = passed / generated
    assert all(math.isfinite(v) for v in values.values())
    rows.append(dict(label=method, display_label=label, n_molecules=23, **values))
with (DATA / 'druglike_means.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
ranks = {k: sorted({r[k] for r in rows}, reverse=k not in LOWER_IS_BETTER)[:2] for k in COLUMNS}

def display(row, key, tex=False):
    value = row[key]
    shown = f'{value:.3f}' if key in LOWER_IS_BETTER else f'{value if key == "num_gen_confs" else value * 100:.1f}'
    # Rank distinct unrounded values; ties share emphasis.
    rank = next((i for i, best in enumerate(ranks[key]) if math.isclose(value, best, rel_tol=0, abs_tol=1e-12)), None)
    if rank is None:
        return shown
    if tex:
        return ('\\textbf{' if rank == 0 else '\\underline{') + shown + '}'
    return '**' + shown + '**' if rank == 0 else '<u>' + shown + '</u>'

note = ('COV-R and COV-P use a strict RMSD < 0.75 Å cutoff; MAT-R and MAT-P are continuous distances. '
        'Coverage, precision, matching distances, and supplied counts are means over 23 molecules, with '
        '2,450 reference records in total. PB pass is the separately measured pooled pass fraction. '
        'Bold and underlining mark the best and second-best distinct unrounded values across the eight methods, '
        'including ties. Higher COV-R, COV-P, and PB pass and lower MAT-R and MAT-P are preferred; supplied '
        'counts are highlighted by size only and describe budget, not quality. '
        'The coverage and distance columns use supplied pools without a common PB filter, and the evaluators '
        'differ in RMSD-failure handling. They remain descriptive, pending harmonized evaluation; unmatched '
        'conformers are not necessarily physically inaccessible. '
        '† Torsion raw starts from the first experimental reference for each molecule and is therefore a '
        'reference-seeded baseline.')
tex_rows = []
md_rows = []
for row in rows:
    seeded = row['label'] == 'torsion_raw_druglike'
    tex_name = row['display_label'] + ('$^{\\dagger}$' if seeded else '')
    tex_rows.append(' & '.join([tex_name] + [display(row, k, True) for k in COLUMNS]) + r' \\')
    md_rows.append('| ' + ' | '.join([row['display_label'] + ('†' if seeded else '')] + [display(row, k) for k in COLUMNS]) + ' |')
tex_note = note.replace('Å', r'\AA{}').replace('†', r'$\dagger$').replace('<', r'$<$')
tex = r'''% Source: results-data/druglike_selected.csv; rendered by build_druglike_table.py.
\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{4pt}
\renewcommand{\arraystretch}{1.15}
\caption{Exploratory coverage and precision for the 23-molecule PLINDER panel.}
\label{tab:druglike}
\par\smallskip
\begin{tabular}{@{}lrrrrrr@{}}
\toprule
Method & \shortstack{Mean\\supplied $N$} & \shortstack{COV-R\\(\%) $\uparrow$} & \shortstack{COV-P\\(\%) $\uparrow$} & \shortstack{MAT-R\\(\AA{}) $\downarrow$} & \shortstack{MAT-P\\(\AA{}) $\downarrow$} & \shortstack{PB pass\\(\%) $\uparrow$} \\
\midrule
''' + '\n'.join(tex_rows) + '\n' + r'''\bottomrule
\end{tabular}
\par\smallskip
\begin{minipage}{\linewidth}
''' + tex_note + '\n' + r'''\end{minipage}
\end{table}
'''
(PAPER / 'tables/druglike.tex').write_text(tex)
md = ('**Table 4. Exploratory coverage and precision for the 23-molecule PLINDER panel.**\n\n'
      '| Method | Mean supplied N | COV-R (%) ↑ | COV-P (%) ↑ | MAT-R (Å) ↓ | MAT-P (Å) ↓ | PB pass (%) ↑ |\n'
      '| --- | --- | --- | --- | --- | --- | --- |\n' + '\n'.join(md_rows) + '\n\n' + note + '\n\n'
      'Sources: [molecule-level records](results_tables/druglike_selected.csv) and '
      '[unrounded table values](results_tables/druglike_means.csv).')
draft = DATA / 'results-working.md'
s, count = re.subn(r'<!-- results-table-5:start -->.*?<!-- results-table-5:end -->',
                  lambda _: '<!-- results-table-5:start -->\n' + md + '\n<!-- results-table-5:end -->',
                  draft.read_text(), flags=re.S)
assert count == 1
draft.write_text(s)
print('Rendered eight Table 4 rows with top-two highlighting.')
