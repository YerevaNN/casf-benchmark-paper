"""Render the landscape figure at its final 6.5-inch manuscript width."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap, to_rgb

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
QWEN = 'qwen_1p7b_v2_coordtoken_step41054'
DATA = json.loads((HERE / 'plot-data.json').read_text())
entries = json.loads((PAPER / 'plinder22-data/entries.json').read_text())
lookup = {(e['method'], e['smiles']): e['key'] for e in entries}
original = np.load(PAPER / 'plinder22-data/nearest-distances.npz', allow_pickle=False)
palette = json.loads((PAPER / 'sampling-data/palette.json').read_text())
colors = {v[0]: v[1] for v in palette.values()}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'text.color': '#29323C', 'axes.labelcolor': '#29323C'})


def render():
    land = json.loads((HERE / 'landscape-data.json').read_text())
    selected = [('loqi', 'loqi_druglike'), ('flowr', 'flowr_druglike'), ('qwen_v2_c17', QWEN)]
    norm = Normalize(0, 3)
    fig = plt.figure(figsize=(6.5, 6.3), layout='constrained')
    grid = fig.add_gridspec(5, 3, height_ratios=[.14, 1, .14, 1, .09])
    for row, name in enumerate(['Staurosporine', 'Actinonin']):
        m = next(m for m in land['molecules'] if m['name'] == name)
        s, d = m['smiles'], m['dedup']
        ix = {k: [i for i, v in enumerate(m['pts']['s']) if m['sources'][v] == k] for k in m['sources']}
        ri = [ix['gt'][j] for j in d['ref_idx']]
        x, y = np.array(d['pts']['x']), np.array(d['pts']['y'])
        all_indices = ri + [i for k, _ in selected for i in ix[k]]
        extent = max(np.ptp(x[all_indices]), np.ptp(y[all_indices])) * 1.15
        cx = (x[all_indices].min() + x[all_indices].max()) / 2
        cy = (y[all_indices].min() + y[all_indices].max()) / 2
        counts = [r['n_clusters'] for r in DATA['clusters'] if r['name'] == name and r['permutation'] > 0]
        heading = fig.add_subplot(grid[row * 2, :])
        heading.axis('off')
        heading.text(0, .5, f'{name} · {len(ri)} references · {np.mean(counts):.1f} clusters at 0.5 Å',
                     va='center', fontsize=9, fontweight='medium')
        for col, (page_key, method) in enumerate(selected):
            ax = fig.add_subplot(grid[row * 2 + 1, col])
            label = DATA['methods'][method]
            rgb = np.array(to_rgb(colors[label]))
            cmap = LinearSegmentedColormap.from_list(method, [.18 * rgb + .82, rgb])
            idx, dist = ix[page_key], d['nn_ref'][page_key]
            assert len(idx) == len(dist)
            r = original[lookup[method, s] + '_reference']
            g = original[lookup[method, s] + '_generated']
            ax.scatter(x[idx], y[idx], c=dist, cmap=cmap, norm=norm,
                       s=3.5, linewidths=0, rasterized=True)
            ax.scatter(x[ri], y[ri], marker='D', c='#333A43', s=12,
                       edgecolors='white', linewidths=.35, zorder=5)
            ax.set(xlim=(cx - extent / 2, cx + extent / 2),
                   ylim=(cy - extent / 2, cy + extent / 2),
                   xticks=[], yticks=[], aspect='equal')
            for spine in ax.spines.values():
                spine.set_color('#BFC6CD')
                spine.set_linewidth(.6)
            ax.set_title(f"{'ABCDEF'[row * 3 + col]}  {label}",
                         loc='left', fontsize=8, pad=5)
            ax.text(.03, .97, f'COV-R {(r <= .75).mean():.1%}\nCOV-P {(g <= .75).mean():.1%}',
                    transform=ax.transAxes, va='top', fontsize=7.5, linespacing=1.3,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .9, 'pad': 1.5})
            ax.text(.03, .03, f'N = {len(idx):,}', transform=ax.transAxes, fontsize=7.5,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .9, 'pad': 1})
            if row == 1:
                cax = fig.add_subplot(grid[4, col])
                bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                                   cax=cax, orientation='horizontal', extend='max')
                bar.set_ticks([0, .75, 1.5, 3])
                bar.ax.tick_params(labelsize=7, length=2)
                bar.set_label('Nearest-reference RMSD (Å)', fontsize=7.5)
    fig.suptitle('Grey diamonds: experimental references', fontsize=8)
    target = PAPER / 'figures/figure-conformer-landscapes'
    for suffix in ['pdf', 'png', 'svg']:
        fig.savefig(target.with_suffix('.' + suffix), dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    render()
