"""Preview aggregate size effects above representation effects; no paper edits."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
summary = pd.read_csv(HERE / 'summary.csv')
selected = summary[summary.tier.isin(['fixed', 'supplied'])]
colors = {'Binned': '#AA589B', 'FSQ': '#0072B2'}
ink, muted = '#29323C', '#626C78'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10,
    'text.color': ink, 'axes.labelcolor': ink,
    'axes.edgecolor': '#89929D', 'axes.linewidth': .7,
    'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.color': muted, 'ytick.color': muted,
    'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'svg.hashsalt': 'qwen-main-two-row',
})
metrics = [
    ('ref', 'best_rmsd', 'CASF ref · 1,044 ligands', 'Mean best RMSD (Å) ↓', (.31, .39)),
    ('plinder', 'best_rmsd', 'PLINDER · 22 molecules', 'Mean MAT-R (Å) ↓', (.37, .45)),
    ('plinder', 'pb', 'PLINDER · 22 molecules', 'PoseBusters pass (%) ↑', (50, 100)),
]
fig, axes = plt.subplots(2, 3, figsize=(12.8, 9.3), sharey='col')
fig.subplots_adjust(left=.085, right=.97, top=.80, bottom=.15, hspace=.90, wspace=.40)
fig.text(.055, .975, 'Model size and representation affect different aspects of generation',
         fontsize=17, va='top')
fig.text(.055, .935, 'The same four configurations, compared by size above and by representation below.',
         fontsize=10, color=muted, va='top')
fig.text(.055, .88, 'Model size: 0.6B → 1.7B', fontsize=14)
fig.text(.055, .455, 'Representation: binned → FSQ', fontsize=14)

rows = []
for col, (cohort, metric, title, ylabel, limits) in enumerate(metrics):
    data = selected[selected.cohort.eq(cohort)].set_index(['representation', 'size'])
    assert len(data) == 4 and data.index.is_unique
    assert np.isfinite(data[metric]).all()
    for row in range(2):
        ax = axes[row, col]
        ax.set_axisbelow(True)
        ax.grid(axis='y', color='#E7EBEE', linewidth=.7)
        ax.set(xlim=(-.25, 1.25), ylim=limits, ylabel=ylabel, xticks=[0, 1],
               xticklabels=['0.6B', '1.7B'] if row == 0 else ['Binned', 'FSQ'])
        ax.set_title(f'{"ABCDEF"[row * 3 + col]}   {title}', loc='left', fontsize=11, pad=12)
        if row == 0:
            for rep in colors:
                ax.plot([0, 1], [data.loc[(rep, size), metric] for size in [.6, 1.7]],
                        color=colors[rep], lw=2.3, zorder=2)
        else:
            # Keep representation colors attached to the endpoints in both rows.
            # Line style and marker distinguish size in this row.
            for size in [.6, 1.7]:
                ax.plot([0, 1], [data.loc[(rep, size), metric] for rep in colors],
                        color='#7C8590', lw=1.7, ls='-' if size == .6 else (0, (4, 3)), zorder=2)
        for rep in colors:
            for size in [.6, 1.7]:
                x = ([.6, 1.7].index(size) if row == 0 else list(colors).index(rep))
                val = float(data.loc[(rep, size), metric])
                other = (data.loc[('FSQ' if rep == 'Binned' else 'Binned', size), metric]
                         if row == 0 else data.loc[(rep, 1.7 if size == .6 else .6), metric])
                assert limits[0] < val < limits[1]
                ax.scatter(x, val, s=77, marker='o' if size == .6 else 'D',
                           facecolors='white' if size == .6 else colors[rep],
                           edgecolors=colors[rep], lw=1.7, zorder=4)
                ax.annotate(f'{val:.1f}%' if metric == 'pb' else f'{val:.3f}',
                            (x, val), xytext=(0, 12 if val > other else -19),
                            textcoords='offset points', ha='center', color=colors[rep], fontsize=10)
                rows.append(dict(panel='ABCDEF'[row * 3 + col], cohort=cohort,
                                 metric=metric, representation=rep, size=size, value=val))

handles = [Line2D([], [], color=color, marker='o', ls='', label=rep)
           for rep, color in colors.items()]
handles += [Line2D([], [], color=muted, marker='o', mfc='white', ls='-', label='0.6B'),
            Line2D([], [], color=muted, marker='D', ls='--', label='1.7B')]
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.52, .073),
           ncol=4, frameon=False, fontsize=10)
fig.text(.055, .043, 'CASF: PB-filtered ensembles, 1,000-candidate target. '
         'PLINDER: full supplied pools; PB measured separately.', fontsize=8, color=muted)
fig.text(.055, .022, 'Matching distances are means across molecules; PB is the pooled pass rate. '
         'Each column uses the same truncated y-axis in both rows.', fontsize=8, color=muted)

stem = HERE / '08-main-size-and-representation'
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(stem.with_suffix('.' + ext), dpi=200, facecolor='white')
plt.close(fig)
pd.DataFrame(rows).to_csv(HERE / f'{stem.name}-values.csv', index=False)
print(f'Saved {stem}.png, .pdf and .svg; checked all 24 plotted values.')
