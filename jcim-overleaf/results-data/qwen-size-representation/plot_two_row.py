"""Preview paired representation and size effects; no manuscript changes."""
from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--simplified', action='store_true',
                    help='Keep marginal percentages and Actinonin; add CASF zooms.')
args = parser.parse_args()
casf = pd.read_csv(HERE / 'casf-per-entry.csv')
plinder = pd.read_csv(HERE / 'plinder-per-molecule.csv')
models = ['binned06', 'binned17', 'fsq06', 'fsq17']
panels = {
    'ref': casf[casf.cohort.eq('ref') & casf.tier.eq('fixed')].pivot(
        index='mol_id', columns='model', values='casf_best_rmsd')[models],
    'plinder': plinder.pivot(index='smiles', columns='model', values='mat_r')[models],
}
names = plinder.drop_duplicates('smiles').set_index('smiles')['name']
assert panels['ref'].shape == (1044, 4)
assert panels['plinder'].shape == (22, 4)
assert all(np.isfinite(w.to_numpy()).all() for w in panels.values())

blue, orange, gray = '#0072B2', '#CE6428', '#9CA5AE'
ink, muted = '#29323C', '#626C78'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10, 'text.color': ink,
    'axes.labelcolor': ink, 'axes.edgecolor': '#89929D',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': .7, 'xtick.color': muted, 'ytick.color': muted,
    'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'svg.hashsalt': 'qwen-representation-size-preview',
})
fig, axes = plt.subplots(2, 2, figsize=(12.4, 12.4))
fig.subplots_adjust(left=.10, right=.97, top=.85, bottom=.14,
                    hspace=.70, wspace=.34)
fig.text(.055, .985, 'How representation and model size change conformer matching',
         fontsize=18, va='top')
fig.text(.055, .951, 'One point per molecule; positive values mean a lower matching distance.',
         fontsize=11, color=muted, va='top')
fig.text(.055, .907, 'Representation: binned → FSQ', fontsize=14, weight='medium')
fig.text(.055, .467, 'Model size: 0.6B → 1.7B', fontsize=14, weight='medium')

statistics, points = [], []
for col, (cohort, w) in enumerate(panels.items()):
    # Preserve the same limits between rows, retaining every outlier.
    limits = (-.78, 1.52) if cohort == 'ref' else (-.30, .68)
    title = 'CASF ref · 1,044 ligands' if cohort == 'ref' else 'PLINDER · 22 molecules'
    effects = [
        ('representation', w.binned06 - w.fsq06, w.binned17 - w.fsq17,
         'FSQ improvement at 0.6B (Å)', 'FSQ improvement at 1.7B (Å)',
         'FSQ closer at both sizes', 'Binned closer at both sizes', '0.6B', '1.7B'),
        ('size', w.binned06 - w.binned17, w.fsq06 - w.fsq17,
         'Size improvement with binned coordinates (Å)', 'Size improvement with FSQ (Å)',
         '1.7B closer with both', '0.6B closer with both', 'Binned', 'FSQ'),
    ]
    for row, (effect, x, y, xlabel, ylabel, good_label, bad_label, xname, yname) in enumerate(effects):
        ax = axes[row, col]
        both = (x > 0) & (y > 0)
        neither = (x < 0) & (y < 0)
        assert min(x.min(), y.min()) > limits[0]
        assert max(x.max(), y.max()) < limits[1]
        colors = np.where(both, blue, np.where(neither, orange, gray))
        ax.set_axisbelow(True)
        ax.grid(color='#E7EBEE', linewidth=.6)
        ax.scatter(x, y, c=colors, s=14 if cohort == 'ref' else 58,
                   alpha=.45 if cohort == 'ref' else .95, edgecolors='none', zorder=3)
        ax.axvline(0, color='#89929D', linewidth=.9, linestyle=(0, (3, 3)))
        ax.axhline(0, color='#89929D', linewidth=.9, linestyle=(0, (3, 3)))
        ax.set(xlim=limits, ylim=limits, xlabel=xlabel, ylabel=ylabel)
        ax.set_aspect('equal', adjustable='box')
        ax.set_title(f'{"ABCD"[row * 2 + col]}   {title}', loc='left', fontsize=12, pad=13)
        if not args.simplified:
            ax.text(.03 if cohort == 'plinder' else .97, .97,
                    f'{good_label}\n{both.sum()}/{len(w)} ({100 * both.mean():.1f}%)',
                    transform=ax.transAxes, ha='left' if cohort == 'plinder' else 'right',
                    va='top', color=blue, fontsize=10)
            ax.text(.03, .03, f'{bad_label}\n{neither.sum()}/{len(w)} ({100 * neither.mean():.1f}%)',
                    transform=ax.transAxes, va='bottom', color=orange, fontsize=9)
        ax.text(.5, -.25,
                f'Improved: {xname} {100 * x.gt(0).mean():.1f}%  ·  {yname} {100 * y.gt(0).mean():.1f}%',
                transform=ax.transAxes, ha='center', color=muted, fontsize=9)
        if args.simplified and cohort == 'ref':
            # Place zooms outside the full panels so no outliers are obscured.
            zoom = ax.inset_axes([1.16, .57, .42, .42])
            zoom.set_axisbelow(True)
            zoom.grid(color='#E7EBEE', linewidth=.5)
            zoom.scatter(x, y, c=colors, s=9, alpha=.55, edgecolors='none', zorder=3)
            for draw in (zoom.axhline, zoom.axvline):
                draw(0, color='#89929D', linewidth=.7, linestyle=(0, (3, 3)))
            zoom.set(xlim=(-.1, .1), ylim=(-.1, .1),
                     xticks=[-.1, 0, .1], yticks=[-.1, 0, .1])
            zoom.set_aspect('equal', adjustable='box')
            zoom.tick_params(labelsize=7, pad=2)
            zoom.set_title('Near zero (±0.10 Å)', fontsize=8, pad=7)
            ax.add_patch(Rectangle((-.1, -.1), .2, .2, fill=False,
                                   edgecolor=ink, linewidth=.8, linestyle=':', zorder=4))
        if cohort == 'plinder':
            offsets = {
                'Actinonin': (-10, -18) if row == 0 else (8, 9),
                'Amprenavir': (8, 7) if row == 0 else (30, 6),
                'Imatinib': (8, -16) if row == 0 else (8, 4),
            }
            if args.simplified:
                offsets = {'Actinonin': offsets['Actinonin']}
            for name, offset in offsets.items():
                identity = names[names.eq(name)].index[0]
                ax.annotate(name, (x[identity], y[identity]), xytext=offset,
                            textcoords='offset points', fontsize=9,
                            ha='right' if offset[0] < 0 else 'left',
                            arrowprops={'arrowstyle': '-', 'lw': .6, 'color': muted}, zorder=4)
        statistics.append(dict(cohort=cohort, effect=effect, n=len(w),
                               improved_x=int(x.gt(0).sum()), improved_y=int(y.gt(0).sum()),
                               improved_both=int(both.sum()), worsened_both=int(neither.sum()),
                               mean_x=float(x.mean()), mean_y=float(y.mean())))
        for identity in w.index:
            points.append(dict(cohort=cohort, effect=effect, identity=identity,
                               name=names[identity] if cohort == 'plinder' else identity,
                               x=x[identity], y=y[identity]))

handles = [Line2D([], [], marker='o', linestyle='', color=color, label=label)
           for color, label in [(blue, 'Improved in both comparisons'),
                                (orange, 'Worsened in both'), (gray, 'Mixed or tied')]]
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.54, .038),
           ncol=3, frameon=False, fontsize=9)
fig.text(.055, .023, 'CASF: per-ligand minimum RMSD after PB filtering, 1,000-candidate target. '
         'PLINDER: per-molecule MAT-R on full supplied pools.', fontsize=8, color=muted)
fig.text(.055, .008, 'Axis ranges are identical between rows within each dataset. '
         + ('Dotted boxes mark the CASF zooms; full panels retain all molecules.' if args.simplified
            else 'All molecules retained; PoseBusters validity is not shown here.'), fontsize=8, color=muted)

stem = HERE / ('07-representation-and-size-simplified' if args.simplified else '06-representation-and-size')
for extension in ['png', 'pdf', 'svg']:
    fig.savefig(stem.with_suffix('.' + extension), dpi=200, facecolor='white')
plt.close(fig)
pd.DataFrame(points).to_csv(HERE / f'{stem.name}-points.csv', index=False)
(HERE / f'{stem.name}-stats.json').write_text(json.dumps(statistics, indent=2) + '\n')
print(pd.DataFrame(statistics).to_string(index=False))
print(f'Saved preview: {stem}.png')
