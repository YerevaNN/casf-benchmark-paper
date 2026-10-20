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
    return method + (r'$^{\dagger}$' if plinder and method == 'Torsion raw' else '')


def write(name, caption, ref, columns, head, body, notes):
    tex = [r'% Rebuild with results-data/three-dataset-tables/build_tables.py.',
           r'\begin{table*}[htbp]', r'\centering\footnotesize',
           r'\setlength{\tabcolsep}{3pt}', r'\renewcommand{\arraystretch}{1.15}',
           r'\caption{' + caption + '}', r'\label{' + ref + '}',
           r'\par\smallskip', r'\resizebox{\linewidth}{!}{%',
           r'\begin{tabular}{' + columns + '}', r'\toprule', *head,
           r'\midrule', *body, r'\bottomrule', r'\end{tabular}%', '}',
           r'\par\smallskip', r'\begin{minipage}{\linewidth}', notes,
           r'\end{minipage}', r'\end{table*}']
    (TABLES / name).write_text('\n'.join(line.rstrip() for line in tex) + '\n')


rank_note = ('Bold and underlining mark the best and second-best distinct unrounded values '
             'among generated methods within each dataset' )
status_note = ('Qwen results precede retraining with the updated benchmark exclusions. ')
plinder_status = ('PLINDER uses the cleaned 22-molecule reference collection. ')
seed_note = (r'$\dagger$ PLINDER torsion sampling starts from an experimental reference. ')
group_head = [
    r'Method & \multicolumn{3}{c}{PLINDER (22)} & \multicolumn{3}{c}{CASF core (94)} & \multicolumn{3}{c}{CASF ref (1,236)} \\',
    r'\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(l){8-10}',
]
columns = r'@{}lrrr@{\hspace{10pt}}rrr@{\hspace{10pt}}rrr@{}'
cohorts = ('plinder', 'core', 'ref')

data = read(1)
methods = list(dict.fromkeys(r['method'] for r in data))
head = group_head + [' & ' + ' & '.join([
    r'\shortstack{Clusters\\0.5 \AA{} $\uparrow$}',
    r'\shortstack{Clusters\\1.0 \AA{} $\uparrow$}',
    r'\shortstack{Largest cluster\\(\%) $\downarrow$}'] * 3) + r' \\']
body = []
for method in methods:
    if method == BASELINE:
        body.append(r'\hdashline')
    values = []
    for cohort in cohorts:
        peers = [r for r in data if r['cohort'] == cohort]
        row = next(r for r in peers if r['method'] == method)
        values += [cell(row, m, peers, pending=cohort == 'plinder' and
                        m == 'largest_percent' and int(row['n_defined']) > 0)
                   for m in ('clusters_0p5', 'clusters_1p0', 'largest_percent')]
    body.append(' & '.join([label(method, True)] + values) + r' \\')
notes = (r'Entries are molecule-level means. Largest cluster denotes the fraction of conformers '
         r'in the largest cluster at 1.0 \AA{}, calculated per molecule before averaging. '
         'PLINDER clustering uses the full supplied pool of conformers per molecule without a common '
         'PoseBusters filter; CASF uses the retained ChEMBL-count ensembles. Absolute cluster counts '
         'therefore reflect different sampling and filtering conditions. CASF means use 92--94 '
         'defined entries on core and 1,220--1,236 on ref, depending on method; PLINDER uses all '
         '22 molecules for each reported method. Retained counts are archived with the source data. '
         'Dashes denote no matched result. ' + rank_note +
         '; ChEMBL3D-PB is unranked. Cluster counts describe geometric diversity, not recovery. ' +
         seed_note + status_note + plinder_status)
write('diversity.tex', 'Geometric diversity across three evaluation sets.', 'tab:diversity', columns, head, body, notes)

data = read(2)
methods = list(dict.fromkeys(r['method'] for r in data))
head = group_head + [
    r' & \shortstack{COV-R@0.5\\(\%) $\uparrow$} & \shortstack{COV-R@0.75\\(\%) $\uparrow$} & \shortstack{MAT-R\\(\AA{}) $\downarrow$} & '
    + ' & '.join([r'\shortstack{Hit@0.5\\(\%) $\uparrow$}', r'\shortstack{Hit@0.75\\(\%) $\uparrow$}',
                  r'\shortstack{Best RMSD\\(\AA{}) $\downarrow$}'] * 2) + r' \\']
body = []
for method in methods:
    if method == BASELINE:
        body.append(r'\hdashline')
    values = []
    for cohort in cohorts:
        peers = [r for r in data if r['cohort'] == cohort]
        row = next(r for r in peers if r['method'] == method)
        values += [cell(row, m, peers, digits=3 if m == 'distance' else 1)
                   for m in ('recovery_0p5', 'recovery_0p75', 'distance')]
    body.append(' & '.join([label(method, True)] + values) + r' \\')
notes = (r'Cutoffs and distances are in \AA{}. PLINDER COV-R is the fraction of experimental references matched at RMSD $\leq t$; '
         'MAT-R averages the nearest-generated RMSD over references. Both are calculated within '
         r'each molecule and then averaged equally over 22 molecules. CASF Hit uses RMSD $\leq t$ '
         'and all 94 core or 1,236 ref entries, counting missing outputs as failures. Best RMSD '
         'averages per-entry minima over defined entries (92--94 on core; 1,220--1,236 on ref). '
         'PLINDER uses supplied pools without a common PoseBusters filter; CASF uses retained '
         'ChEMBL-count ensembles. Dashes denote no matched result. PLINDER values use the same '
         r'cleaned reference distances as Table~\ref{tab:druglike}; all nearest distances are finite. ' +
         rank_note + '; ChEMBL3D-PB is unranked. ' + seed_note + status_note + plinder_status)
write('count-recovery.tex', 'Recovery of observed bound conformations across three evaluation sets.',
      'tab:count-recovery', columns, head, body, notes)

data = read(3)
fixed = [r for r in data if r['cohort'] == 'core' and r['budget'] == 'fixed']
methods = [r['method'] for r in sorted(fixed, key=lambda r: number(r, 'distance'), reverse=True)] + [BASELINE]
head = [r'Method & Candidate target & \multicolumn{3}{c}{CASF core (94)} & \multicolumn{3}{c}{CASF ref (1,236)} \\',
        r'\cmidrule(lr){3-5}\cmidrule(l){6-8}',
        ' & & ' + ' & '.join([r'\shortstack{Best RMSD\\(\AA{}) $\downarrow$}',
                             r'\shortstack{Hit@0.75\\(\%) $\uparrow$}',
                             r'\shortstack{Clusters\\0.5 \AA{} $\uparrow$}'] * 2) + r' \\']
body = []
for method in methods:
    if method == BASELINE:
        body.append(r'\hdashline')
    for i, budget in enumerate(('chembl_count',) if method == BASELINE else ('fixed', 'chembl_count')):
        values = []
        for cohort in ('core', 'ref'):
            peers = [r for r in data if r['cohort'] == cohort and r['budget'] == budget]
            row = next(r for r in peers if r['method'] == method)
            values += [cell(row, m, peers, digits=3 if m == 'distance' else 1)
                       for m in ('distance', 'recovery', 'clusters')]
        target = 'Stored' if method == BASELINE else '1,000' if budget == 'fixed' else 'ChEMBL-count'
        body.append(' & '.join([method if i == 0 else '', target] + values) + r' \\')
    if method != BASELINE:
        body.append(r'\addlinespace[2pt]')
notes = (r'Targets precede PoseBusters filtering. Hit@0.75 uses all 94 core or 1,236 ref entries '
         'and counts missing outputs as failures. Best RMSD and cluster counts use defined entries; '
         'method-specific denominators and retained counts are archived with the source data. '
         'ChEMBL3D-PB is the same stored ensemble and has no 1,000-candidate result. Generators '
         'follow decreasing unrounded core Best RMSD at 1,000 candidates. ' + rank_note +
         ' and candidate target; ChEMBL3D-PB is unranked. Targets do not imply equal retained '
         'counts or computational costs. ' + status_note)
write('recovery.tex', 'Recovery and geometric diversity at two candidate targets on CASF core and ref.',
      'tab:recovery', r'@{}llrrr@{\hspace{12pt}}rrr@{}', head, body, notes)
