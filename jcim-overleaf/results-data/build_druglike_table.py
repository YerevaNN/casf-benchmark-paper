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
    ('rdkit_random_minimized_22_druglike', 'RDKit minimized'),
    ('torsion_raw_druglike', 'Torsion raw'),
    ('torsion_minimized_22_druglike', 'Torsion minimized'),
    ('loqi_druglike', 'LoQI'),
    ('torsional_diffusion_druglike', 'Torsional Diffusion'),
    ('mcf_drugs_l_druglike', 'MCF drugs-L'),
    ('nextmol_dmt_l_druglike', 'NExT-Mol DMT-L'),
    ('flowr_druglike', 'FlowR'),
    ('qwen_1p7b_fsq_bigdata_step47023', 'Qwen 1.7B FSQ'),
]
COLUMNS = ['num_gen_confs', 'cov_r_075', 'cov_p_075', 'mat_r', 'mat_p', 'pooled_pb_pass']
DISPLAY_COLUMNS = ['cov_r_075', 'cov_p_075', 'mat_r', 'mat_p', 'num_gen_confs', 'pooled_pb_pass']
LOWER_IS_BETTER = {'mat_r', 'mat_p'}
with (DATA / 'druglike_selected.csv').open() as handle:
    records = list(csv.DictReader(handle))
assert {r['label'] for r in records} == {m for m, _ in METHODS}
rows = []
identities = None
for method, label in METHODS:
    selected = [r for r in records if r['label'] == method]
    current = {r['smiles'] for r in selected}
    assert len(selected) == len(current) == 22
    if identities is None:
        identities = current
    assert current == identities
    assert sum(int(r['num_true_confs']) for r in selected) == 1067
    values = {k: mean(float(r[k]) for r in selected) for k in COLUMNS[:-1]}
    generated = sum(int(r['n_generated']) for r in selected)
    passed = sum(int(r['n_pb_pass']) for r in selected)
    assert generated > 0 and 0 <= passed <= generated
    values['pooled_pb_pass'] = passed / generated
    assert all(math.isfinite(v) for v in values.values())
    rows.append(dict(label=method, display_label=label, n_molecules=22, **values))
with (DATA / 'druglike_means.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
ranks = {k: sorted({r[k] for r in rows}, reverse=k not in LOWER_IS_BETTER)[:2] for k in COLUMNS}

def display(row, key, tex=False):
    value = row[key]
    shown = f'{value:.3f}' if key in LOWER_IS_BETTER else f'{value if key == "num_gen_confs" else value * 100:.1f}'
    if key == 'num_gen_confs':
        return shown
    # Rank distinct unrounded values; ties share emphasis.
    rank = next((i for i, best in enumerate(ranks[key]) if math.isclose(value, best, rel_tol=0, abs_tol=1e-12)), None)
    if rank is None:
        return shown
    if tex:
        return ('\\textbf{' if rank == 0 else '\\underline{') + shown + '}'
    return '**' + shown + '**' if rank == 0 else '<u>' + shown + '</u>'

note = ('COV-R/P use RMSD ≤ 0.75 Å; MAT-R/P are matching distances. Coverage, precision, distances, '
        'and supplied counts are equal-weight means over 22 molecules with 1,067 cleaned references. '
        'Full supplied pools are used without PB filtering; PB pass is measured separately and pooled '
        'across conformers. Bold/underlining mark the best/second-best distinct unrounded performance '
        'values, sharing ties. Supplied counts are descriptive and unranked. '
        '† Both torsion baselines are seeded from the first experimental reference in the original collection.')
tex_rows = []
md_rows = []
for index, row in enumerate(rows):
    seeded = row['label'] in {'torsion_raw_druglike', 'torsion_minimized_22_druglike'}
    tex_name = row['display_label'] + ('$^{\\dagger}$' if seeded else '')
    if index % 2 == 0:
        tex_rows.append(r'\rowcolor{tableShade}')
    tex_rows.append(' & '.join([tex_name] + [display(row, k, True) for k in DISPLAY_COLUMNS]) + r' \\')
    md_rows.append('| ' + ' | '.join([row['display_label'] + ('†' if seeded else '')] + [display(row, k) for k in DISPLAY_COLUMNS]) + ' |')
tex_note = note.replace('Å', r'\AA{}').replace('†', r'$\dagger$').replace('≤', r'$\leq$')
tex = r'''% Source: results-data/druglike_selected.csv; rendered by build_druglike_table.py.
\begin{table*}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{2pt}
\renewcommand{\arraystretch}{1.35}
\caption{Exploratory coverage and precision for the 22-molecule PLINDER panel.}
\label{tab:druglike}
\par\smallskip
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrrrr@{\hspace{10pt}}rr@{}}
\toprule
 & \multicolumn{4}{c}{\textcolor{plinderTeal}{\textbf{Reference agreement}}} & \multicolumn{2}{c}{\textbf{Pool checks}} \\
\cmidrule(lr){2-5}\cmidrule(l){6-7}
\textbf{Method} & \shortstack{COV-R\\(\%) $\uparrow$} & \shortstack{COV-P\\(\%) $\uparrow$} & \shortstack{MAT-R\\(\AA{}) $\downarrow$} & \shortstack{MAT-P\\(\AA{}) $\downarrow$} & Supplied $N$ & \shortstack{PB pass\\(\%) $\uparrow$} \\
\midrule
''' + '\n'.join(tex_rows) + '\n' + r'''\bottomrule
\end{tabular*}
\par\smallskip
\begin{minipage}{\linewidth}
''' + tex_note + '\n' + r'''\end{minipage}
\end{table*}
'''
(PAPER / 'tables/druglike.tex').write_text(tex)
# Supporting Table S4 uses the same ten methods and unrounded means.
support = r"""\begin{table}[htbp]
\centering\small
\caption{Matching distances for the cleaned 22-molecule PLINDER panel.
MAT-R averages nearest generated-conformer distances over references; MAT-P averages
nearest reference distances over generated conformers. Values are equal-weight molecule means.
Supplied pools are not filtered by PoseBusters; all nearest distances are finite.
Both torsion baselines are seeded from an experimental reference.}
\label{tab:supporting-druglike}
\begin{tabular}{lrr}
\toprule
Method & MAT-R (\AA{}) & MAT-P (\AA{}) \\
\midrule
"""
for row in rows:
    name = row['display_label'] + (' (reference-seeded)' if row['label'] in {'torsion_raw_druglike', 'torsion_minimized_22_druglike'} else '')
    support += f"{name} & {row['mat_r']:.3f} & {row['mat_p']:.3f}" + r' \\' + '\n'
support += r"""\bottomrule
\end{tabular}
\end{table}
"""
(PAPER / 'tables/supporting-druglike.tex').write_text(support)
print('Rendered Tables 4 and S4 from the unified PLINDER-22 source.')
