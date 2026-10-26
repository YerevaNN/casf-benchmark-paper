"""Render Supporting Table S9 from the common-protocol energy summaries."""
from pathlib import Path
import math
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
core = pd.read_csv(PAPER / 'energy-data/summary.csv')
ref = pd.read_csv(HERE / 'ref/summary.csv').set_index('method').loc[core.method].reset_index()


def number(value, signed=False):
    # Keep extreme finite values visible without widening the table excessively.
    if abs(value) >= 1e6:
        exponent = math.floor(math.log10(abs(value)))
        return rf'${value / 10**exponent:.1f}\times10^{{{exponent}}}$'
    return f'{value:+.1f}' if signed else f'{value:.1f}'


lines = [r'''\begin{table}[p]
\centering\footnotesize
\setlength{\tabcolsep}{3pt}
\renewcommand{\arraystretch}{1.12}
\caption{Energy summaries after common hydrogen preparation at the ChEMBL-count target, separately for CASF core and ref. All energy values are in kcal/mol.}
\label{tab:energy-summary}
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrrrrr@{}}
\toprule
Method & $n$ & \shortstack{Median\\energy SD} & \shortstack{Mean\\energy SD} & \shortstack{Median\\mean energy} & \shortstack{Median paired\\energy difference} \\
\midrule''']
for label, data in [('CASF core', core), ('CASF ref', ref)]:
    lines.append(r'\multicolumn{6}{l}{\textbf{' + label + r'}} \\')
    lines.append(r'\midrule')
    for row in data.itertuples():
        if row.method == 'chembl3d_gt_pb':
            lines.append(r'\hdashline')
        values = [row.label, str(row.n_molecules), number(row.median_energy_std),
                  number(row.mean_energy_std), number(row.median_mean_energy),
                  number(row.median_mean_energy_delta_chembl, signed=True)]
        lines.append(' & '.join(values) + r' \\')
    lines.append(r'\midrule' if label == 'CASF core' else r'\bottomrule')
lines.append(r'''\end{tabular*}
\par\smallskip
\begin{minipage}{\linewidth}
$n$: ensembles included. Energy SD: within-ensemble population standard deviation. Paired energy differences are relative to ChEMBL3D-PB for the same molecule; positive values indicate higher mean energy.

Finite energies were obtained for all 77,986 retained core conformers and 835,972 of 835,973 retained ref conformers. The undefined energy belongs to FlowR. The large ref LoQI mean SD includes finite extreme values.
\end{minipage}
\end{table}''')
(PAPER / 'tables/energy-summary.tex').write_text('\n'.join(lines) + '\n')
