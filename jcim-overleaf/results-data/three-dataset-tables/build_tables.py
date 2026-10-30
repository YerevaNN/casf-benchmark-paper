"""Render Tables 1–3 from the archived CASF and refreshed PLINDER-22 values."""
import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TABLES = HERE.parents[1] / 'tables'
BASELINE = 'ChEMBL3D-PB'


def read(number):
    with (HERE / f'table{number}.csv').open() as stream:
        return list(csv.DictReader(stream))


def number(row, metric):
    return float(row[metric]) if row.get(metric) else math.nan


def cell(row, metric, peers, digits=1, pending=False):
    value = number(row, metric)
    if math.isnan(value):
        return r'\textit{pending}' if pending else '---'
    text = f'{value:.{digits}f}'
    if row['method'] == BASELINE:
        return text
    lower = metric in ('distance', 'largest_percent')
    ranks = sorted({number(p, metric) for p in peers
                    if p['method'] != BASELINE and math.isfinite(number(p, metric))},
                   reverse=not lower)[:2]
    if value == ranks[0]:
        return r'\textbf{' + text + '}'
    if len(ranks) > 1 and value == ranks[1]:
        return r'\underline{' + text + '}'
    return text


def label(method, plinder):
    return method + (r'$^{\dagger}$' if plinder and method in {'Torsion raw', 'Torsion minimized'} else '')


def write(name, caption, ref, columns, head, body, notes):
    tex = [r'% Rebuild with results-data/three-dataset-tables/build_tables.py.',
           r'\begin{table*}[htbp]', r'\centering\footnotesize',
           r'\setlength{\tabcolsep}{2pt}', r'\renewcommand{\arraystretch}{1.35}',
           r'\caption{' + caption + '}', r'\label{' + ref + '}',
           r'\par\smallskip',
           r'\begin{tabular*}{\linewidth}{' + columns + '}', r'\toprule', *head,
           r'\midrule', *body, r'\bottomrule', r'\end{tabular*}',
           r'\par\smallskip', r'\begin{minipage}{\linewidth}', notes,
           r'\end{minipage}', r'\end{table*}']
    (TABLES / name).write_text('\n'.join(line.rstrip() for line in tex) + '\n')


rank_note = 'Bold/underlined: best/second best. '
seed_note = r'$\dagger$ Reference-seeded. '
group_head = [
    r' & \multicolumn{3}{c}{\textcolor{plinderTeal}{\textbf{PLINDER (22)}}} & \multicolumn{3}{c}{\textbf{CASF core (94)}} & \multicolumn{3}{c}{\textbf{CASF ref (1,044)}} \\',
    r'\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(l){8-10}',
]
columns = r'@{\extracolsep{\fill}}lrrr@{\hspace{8pt}}rrr@{\hspace{8pt}}rrr@{}'
cohorts = ('plinder', 'core', 'ref')

data = read(1)
methods = list(dict.fromkeys(r['method'] for r in data))
head = group_head + [r'\textbf{Method} & ' + ' & '.join([
    r'C0.5 $\uparrow$', r'C1.0 $\uparrow$', r'Largest $\downarrow$'] * 3) + r' \\']
body = []
for index, method in enumerate(methods):
    if method == BASELINE:
        body.append(r'\hdashline')
    elif index % 2 == 0:
        body.append(r'\rowcolor{tableShade}')
    values = []
    for cohort in cohorts:
        peers = [r for r in data if r['cohort'] == cohort]
        row = next(r for r in peers if r['method'] == method)
        values += [cell(row, m, peers, pending=cohort == 'plinder' and
                        m == 'largest_percent' and int(row['n_defined']) > 0)
                   for m in ('clusters_0p5', 'clusters_1p0', 'largest_percent')]
    body.append(' & '.join([label(method, True)] + values) + r' \\')
notes = (r'C0.5/C1.0: mean cluster counts at 0.5/1.0 \AA{}. Largest: mean largest-cluster '
         r'share at 1.0 \AA{} (\%). CASF uses the ChEMBL-count target. '
         'Dashes: unavailable. ' + rank_note + seed_note)
write('diversity.tex', 'Geometric diversity across three evaluation sets.', 'tab:diversity', columns, head, body, notes)

data = read(2)
methods = list(dict.fromkeys(r['method'] for r in data))
head = group_head + [
    r'\textbf{Method} & \shortstack{COV-R\\@0.5 $\uparrow$} & \shortstack{COV-R\\@0.75 $\uparrow$} & MAT-R $\downarrow$ & '
    + ' & '.join([r'\shortstack{Hit\\@0.5 $\uparrow$}', r'\shortstack{Hit\\@0.75 $\uparrow$}',
                  r'\shortstack{Best\\RMSD $\downarrow$}'] * 2) + r' \\']
body = []
for index, method in enumerate(methods):
    if method == BASELINE:
        body.append(r'\hdashline')
    elif index % 2 == 0:
        body.append(r'\rowcolor{tableShade}')
    values = []
    for cohort in cohorts:
        peers = [r for r in data if r['cohort'] == cohort]
        row = next(r for r in peers if r['method'] == method)
        values += [cell(row, m, peers, digits=3 if m == 'distance' else 1)
                   for m in ('recovery_0p5', 'recovery_0p75', 'distance')]
    body.append(' & '.join([label(method, True)] + values) + r' \\')
notes = (r'COV-R/Hit are percentages; cutoffs and distances are in \AA{}. '
         'CASF uses the ChEMBL-count target. Dashes: unavailable. ' + rank_note + seed_note)
write('count-recovery.tex', 'Recovery of observed bound conformations across three evaluation sets.',
      'tab:count-recovery', columns, head, body, notes)

data = read(3)
fixed = [r for r in data if r['cohort'] == 'core' and r['budget'] == 'fixed']
methods = [r['method'] for r in sorted(fixed, key=lambda r: number(r, 'distance'), reverse=True)] + [BASELINE]
head = [r' & & \multicolumn{3}{c}{\textbf{CASF core (94)}} & \multicolumn{3}{c}{\textbf{CASF ref (1,044)}} \\',
        r'\cmidrule(lr){3-5}\cmidrule(l){6-8}',
        r'\textbf{Method} & \textbf{Target} & ' + ' & '.join([r'Best RMSD $\downarrow$',
                             r'Hit@0.75 $\uparrow$', r'C0.5 $\uparrow$'] * 2) + r' \\']
body = []
for index, method in enumerate(methods):
    if method == BASELINE:
        body.append(r'\hdashline')
    for i, budget in enumerate(('chembl_count',) if method == BASELINE else ('fixed', 'chembl_count')):
        if method != BASELINE and index % 2 == 0:
            body.append(r'\rowcolor{tableShade}')
        values = []
        for cohort in ('core', 'ref'):
            peers = [r for r in data if r['cohort'] == cohort and r['budget'] == budget]
            row = next(r for r in peers if r['method'] == method)
            values += [cell(row, m, peers, digits=3 if m == 'distance' else 1)
                       for m in ('distance', 'recovery', 'clusters')]
        target = 'Stored' if method == BASELINE else '1,000' if budget == 'fixed' else 'ChEMBL'
        body.append(' & '.join([method if i == 0 else '', target] + values) + r' \\')
    if method != BASELINE:
        body.append(r'\addlinespace[2pt]')
notes = (r'ChEMBL: stored-conformer-count target. Best RMSD is in \AA{}; '
         r'Hit@0.75 is a percentage; C0.5 is the mean cluster count at 0.5 \AA{}. ' + rank_note)
write('recovery.tex', 'Recovery and geometric diversity at two candidate targets on CASF core and ref.',
      'tab:recovery', r'@{\extracolsep{\fill}}llrrr@{\hspace{8pt}}rrr@{}', head, body, notes)
