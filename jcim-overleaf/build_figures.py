"""Recreate manuscript figures from the public dashboard release.

Run from the project environment: python jcim-overleaf/build_figures.py
Use --refresh-data to verify local databases against public-release.json and
re-export their underlying rows. Default rendering uses the archived CSVs only.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sqlite3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'figure-data'
OUT = ROOT / 'figures'
QWEN = 'qwen_1p7b_fsq_bigdata_step47023'
# Shared colors plus distinct shapes: no scientific distinction relies on color alone.
STYLE = {
    QWEN: ('Qwen 1.7B FSQ', '#0072B2', 'o'),
    'loqi_raw': ('LoQI', '#D55E00', 's'),
    'flowr_raw': ('FlowR', '#009E73', '^'),
    'torsional_diffusion_raw': ('Torsional Diffusion', '#CC79A7', 'D'),
    'mcf_drugs_l_raw': ('MCF drugs-L', '#A88700', 'v'),
    'nextmol_dmt_l_raw': ('NExT-Mol DMT-L', '#56B4E9', 'P'),
    'rdkit_random_raw': ('RDKit raw', '#737373', 'X'),
    'torsion_raw': ('Torsion raw', '#8C613C', '*'),
    'chembl3d_gt_pb': ('ChEMBL3D-PB', '#252525', 'h'),
}
KEYS = list(STYLE)
STYLE.update({
    'rdkit_random_minimized': ('RDKit minimized', '#999999', 'd'),
    'torsion_minimized': ('Torsion minimized', '#B39B72', 'p'),
})
PANEL = ['rdkit_random_raw', 'rdkit_random_minimized', 'torsion_raw',
         'torsion_minimized', 'loqi_raw', 'torsional_diffusion_raw',
         'mcf_drugs_l_raw', 'nextmol_dmt_l_raw', 'flowr_raw', QWEN]
GALLERY = None
FOCUS = [QWEN, 'loqi_raw', 'flowr_raw', 'chembl3d_gt_pb']
DRUG = {QWEN: QWEN, 'loqi_raw': 'loqi_druglike', 'flowr_raw': 'flowr_druglike',
        'torsional_diffusion_raw': 'torsional_diffusion_druglike',
        'mcf_drugs_l_raw': 'mcf_drugs_l_druglike', 'nextmol_dmt_l_raw': 'nextmol_dmt_l_druglike'}
THRESHOLDS = [('0p25', .25), ('0p5', .5), ('0p75', .75), ('2p0', 2.)]
RADII = [('0p5', .5), ('1p0', 1.), ('2p0', 2.), ('3p0', 3.)]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
    'axes.labelsize': 8, 'axes.titlesize': 9, 'axes.titleweight': 'normal',
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.linewidth': .65, 'lines.linewidth': 1.3, 'lines.markersize': 4,
    'xtick.major.width': .6, 'ytick.major.width': .6,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'savefig.facecolor': 'white', 'axes.spines.top': False, 'axes.spines.right': False})


def table(db, name):
    with sqlite3.connect(f'file:{db}?mode=ro', uri=True) as con:
        return pd.read_sql_query(f'SELECT * FROM {name}', con)


def refresh_data():
    release = json.loads((DATA / 'public-release.json').read_text())
    verified = []
    for asset in release['assets']:
        path = ROOT.parent / 'data/results' / asset['name']
        digest = 'sha256:' + hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == asset['digest'], f'Public release mismatch: {path}'
        verified.append({'name': path.name, 'sha256': digest, 'url': asset['browser_download_url']})
    main = ROOT.parent / 'data/results/casf_analysis_dashboard.sqlite'
    extended = ROOT.parent / 'data/results/extended_casf_analysis.sqlite'
    frame = table(main, 'per_ligand_long')
    sums, strata, raw = [], [], []
    for cohort, expected_n in [('core', 94), ('ref', 1236)]:
        base = frame[(frame.ligand_set == cohort) & (frame.method == 'chembl3d_gt_pb')].set_index('mol_id')
        assert base.index.is_unique and len(base) == expected_n
        size = pd.cut(base.heavy_atoms, [0, 20, 30, 40, np.inf], right=False, labels=['<20', '20–29', '30–39', '≥40'])
        rotors = pd.cut(base.rotatable_bonds, [-1, 3, 6, 8, np.inf], labels=['0–3', '4–6', '7–8', '≥9'])
        assert size.notna().all() and rotors.notna().all()
        for key in KEYS:
            for tier in (['reference'] if key == 'chembl3d_gt_pb' else ['chembl_count', 'dynamic', 'fixed']):
                method = key if tier == 'reference' else f'{key}_{tier}'
                rows = frame[(frame.ligand_set == cohort) & (frame.method == method)].set_index('mol_id')
                assert rows.index.is_unique and set(rows.index) == set(base.index), (cohort, method)
                rows = rows.reindex(base.index)
                row = dict(cohort=cohort, key=key, method=method, tier=tier, n=expected_n,
                           evaluated=int(rows.casf_best_rmsd.notna().sum()),
                           mean_retained=rows.post_pb_confs.fillna(0).mean())
                for tag, _ in THRESHOLDS:
                    row['hit_' + tag] = rows['casf_hit_' + tag].fillna(0).mean() * 100
                for tag, _ in RADII:
                    row['clusters_' + tag] = rows['greedy_clusters_' + tag].mean()
                sums.append(row)
                if tier not in ('fixed', 'reference'):
                    continue
                cols = ['casf_best_rmsd', 'casf_hit_0p75', 'post_pb_confs']
                entry = rows[cols].copy()
                entry['heavy_atoms_reference'] = base.heavy_atoms
                entry['rotors_reference'] = base.rotatable_bonds
                entry['cohort'], entry['key'] = cohort, key
                raw.append(entry.reset_index())
                for descriptor, bins in [('size', size), ('rotors', rotors)]:
                    for label in bins.cat.categories:
                        selected = rows[bins == label]
                        strata.append(dict(cohort=cohort, key=key, descriptor=descriptor, group=label,
                            n=len(selected), hits=int(selected.casf_hit_0p75.fillna(0).sum()),
                            hit=selected.casf_hit_0p75.fillna(0).mean() * 100))
    summary = pd.DataFrame(sums)
    summary.to_csv(DATA / 'casf_summary.csv', index=False)
    pd.DataFrame(strata).to_csv(DATA / 'casf_strata.csv', index=False)
    pd.concat(raw).to_csv(DATA / 'casf_selected_entries.csv', index=False)
    # Preserve the public displayed aggregates as an audit: they omit undefined hits.
    display = table(main, 'comparison_rows')
    display[display.method.isin(summary.method)].to_csv(DATA / 'dashboard_displayed_aggregates.csv', index=False)
    drug = table(extended, 'extended_druglike_per_molecule')
    drug = drug[drug.label.isin(DRUG.values())].copy()
    assert drug.groupby('label').size().eq(23).all()
    drug[['label', 'smiles', 'name', 'num_true_confs', 'num_gen_confs', 'cov_r_075', 'cov_p_075', 'mat_r', 'mat_p', 'pb_pass_rate']].to_csv(DATA / 'drug_molecule_metrics.csv', index=False)
    energy = table(extended, 'extended_energy_window_summary')
    energy.to_csv(DATA / 'historical_energy_summary.csv', index=False)
    (DATA / 'provenance.json').write_text(json.dumps({
        'release': release['tag_name'], 'release_url': release['html_url'],
        'published_at': release['published_at'], 'verified_on': '2026-09-28',
        'assets': verified, 'qwen_main_checkpoint': QWEN,
        'hit_denominator': 'All mapped entries: core 94; ref 1236; missing outputs are misses.',
        'cluster_denominator': 'Entries with defined clustering results.',
        'budget': 'Candidate-tier endpoints, not random fixed-valid-K curves.',
        'energy_status': 'Historical September 14 sidecar; diagnostic only, excluded from manuscript.',
        'drug_status': 'Supplied pools; no common PB filter; failure policies differ.',
        'method_selection': 'Main CASF figures include all ten selected pipelines plus the stored ChEMBL3D-PB comparison; core only, recovery at 0.75 Å. Additional Qwen checkpoints are outside the main panel.',
    }, indent=2) + '\n')


def style_axis(ax, panel, title, grid='y'):
    ax.set_title(title, loc='left', pad=12)
    ax.text(-.13, 1.055, panel, transform=ax.transAxes, fontsize=11, weight='bold')
    ax.grid(axis=grid, color='#E5E7EB', linewidth=.55)
    ax.set_axisbelow(True)
    ax.tick_params(length=3)


def legend(fig, keys=KEYS, ncol=3, y=.015):
    handles = [Line2D([0], [0], color=STYLE[k][1], marker=STYLE[k][2],
               ls='--' if k == 'chembl3d_gt_pb' else '-', label=STYLE[k][0]) for k in keys]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, y), ncol=ncol,
               frameon=False, handlelength=2.1, columnspacing=1.6)


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    # Fixed 7-inch canvas, embedded vector fonts; 450 dpi color raster companion.
    fig.savefig(OUT / f'{name}.pdf', metadata={'Title': name, 'Creator': 'casf-benchmark manuscript/build_figures.py'})
    svg = OUT / f'{name}.svg'
    fig.savefig(svg)
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
    fig.savefig(OUT / f'{name}.png', dpi=450)
    if GALLERY is not None:
        GALLERY.savefig(fig)
    plt.close(fig)


def get_row(summary, cohort, key, tier='fixed'):
    if key == 'chembl3d_gt_pb':
        tier = 'reference'
    selected = summary[(summary.cohort == cohort) & (summary.key == key) & (summary.tier == tier)]
    assert len(selected) == 1
    return selected.iloc[0]


def recovery(tiers):
    core = tiers[tiers.cohort == 'core']
    fig, axes = plt.subplots(1, 2, figsize=(7, 4.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.25, right=.97, top=.86, bottom=.22, wspace=.20)
    ref = core[core.key == 'chembl3d_gt_pb'].iloc[0]
    for ax, tier, panel, title in zip(axes, ['chembl_count', 'fixed'], ['A', 'B'],
                                    ['ChEMBL-count target', '1,000-candidate target']):
        for y, key in enumerate(PANEL):
            _, color, marker = STYLE[key]
            row = core[(core.key == key) & (core.tier == tier)].iloc[0]
            ax.plot(row.hit_0p75, y, marker=marker, color=color, ls='none', ms=5)
            ax.annotate(f'{row.hit_0p75:.1f}', (row.hit_0p75, y), xytext=(6, 0),
                        textcoords='offset points', va='center', fontsize=7,
                        bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': .3})
        ax.axvline(ref.hit_0p75, color='#252525', ls='--', lw=.9)
        ax.set(xlabel='Recovery at 0.75 Å (%)', xlim=(50, 100),
               ylim=(len(PANEL)-.4, -.7))
        style_axis(ax, panel, title, grid='x')
    axes[0].set_yticks(range(len(PANEL)), [STYLE[k][0] for k in PANEL])
    axes[1].tick_params(axis='y', left=False)
    fig.legend(handles=[Line2D([], [], color='#252525', ls='--',
                              label=f'Stored ChEMBL3D-PB: {ref.hit_0p75:.1f}% recovery')],
               loc='lower center', bbox_to_anchor=(.60, .075), frameon=False, fontsize=7)
    fig.text(.60, .025, 'Core: 94 ligands · candidate targets apply before PoseBusters filtering.',
             ha='center', fontsize=7, color='#555555')
    save(fig, 'core-recovery-budget')


def diversity(comparison):
    fig, ax = plt.subplots(figsize=(7, 3.7))
    fig.subplots_adjust(left=.10, right=.97, top=.87, bottom=.30)
    offsets = {QWEN: (6, 6), 'loqi_raw': (-5, 9), 'flowr_raw': (-28, 9),
               'torsion_raw': (6, -12), 'chembl3d_gt_pb': (6, 6)}
    for key in PANEL + ['chembl3d_gt_pb']:
        row = comparison[comparison.key == key].iloc[0]
        name, color, marker = STYLE[key]
        ax.scatter(row.clusters_1p0, row.hit_0p75, marker=marker, s=36,
                   c=color, edgecolors='white', linewidths=.4, zorder=3)
        if key in offsets:
            ax.annotate(name, (row.clusters_1p0, row.hit_0p75), xytext=offsets[key],
                        textcoords='offset points', fontsize=7, color=color)
    ax.set(xlabel='Mean cluster count at 1.0 Å', ylabel='Recovered ligands at 0.75 Å (%)',
           xlim=(0, 125), ylim=(75, 96))
    style_axis(ax, '', 'Diversity and recovery · core (n = 94)')
    legend(fig, PANEL + ['chembl3d_gt_pb'], ncol=4, y=.01)
    save(fig, 'diversity-recovery')


def size_figure(strata):
    keys = PANEL + ['chembl3d_gt_pb']
    fig, axes = plt.subplots(1, 2, figsize=(7, 5.1), sharey=True)
    fig.subplots_adjust(left=.24, right=.98, top=.84, bottom=.28, wspace=.15)
    for j, (descriptor, groups) in enumerate([
        ('size', ['<20', '20–29', '30–39', '≥40']),
        ('rotors', ['0–3', '4–6', '7–8', '≥9'])]):
        ax = axes[j]
        selected = strata[strata.descriptor == descriptor]
        values = selected.pivot(index='key', columns='group', values='hit').loc[keys, groups]
        hits = selected.pivot(index='key', columns='group', values='hits').loc[keys, groups]
        counts = selected.pivot(index='key', columns='group', values='n').loc[keys, groups]
        assert (counts == counts.iloc[0]).all().all()
        mesh = ax.imshow(values, cmap='Blues', vmin=0, vmax=100, aspect='auto')
        for y, key in enumerate(keys):
            for x, group in enumerate(groups):
                ax.text(x, y, f'{hits.loc[key, group]}/{counts.loc[key, group]}',
                        ha='center', va='center', fontsize=7,
                        color='white' if values.loc[key, group] >= 65 else '#202C34')
        ns = counts.iloc[0].tolist()
        ax.set_xticks(range(4), [f'{group}'+('*' if n < 5 else '') for group, n in zip(groups, ns)])
        ax.set_yticks(range(len(keys)), [STYLE[k][0] for k in keys])
        ax.set_xlabel('Heavy atoms' if descriptor == 'size' else 'Rotatable bonds')
        ax.set_title(('A   Molecular size' if j == 0 else 'B   Flexibility'), loc='left', pad=12)
        ax.tick_params(length=0)
        ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(keys), 1), minor=True)
        ax.grid(which='minor', color='white', linewidth=1)
        ax.tick_params(which='minor', bottom=False, left=False)
        for x, n in enumerate(ns):
            if n < 5:
                ax.add_patch(Rectangle((x-.5, -.5), 1, len(keys), fill=False,
                                       edgecolor='#D55E00', linestyle='--', linewidth=1.2))
        for spine in ax.spines.values():
            spine.set_visible(False)
    bar_ax = fig.add_axes([.40, .14, .43, .025])
    fig.colorbar(mesh, cax=bar_ax, orientation='horizontal', label='Recovered ligands at 0.75 Å (%)')
    fig.text(.60, .95, 'Cells show recovered / evaluated ligands in core (n = 94).', ha='center', fontsize=7)
    fig.text(.60, .015, '* Fewer than five ligands: descriptive counts, not stable rankings.',
             ha='center', fontsize=7, color='#555555')
    save(fig, 'size-flexibility')


def drug_figure(drug):
    rng = np.random.default_rng(20260924)
    keys = [key for key in PANEL if key in DRUG]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7, 4.05), gridspec_kw={'width_ratios':[1,1]})
    fig.subplots_adjust(left=.19, right=.98, top=.85, bottom=.24, wspace=.37)
    exports = []
    # Preserve the archived bootstrap draw order when changing display order.
    for key in DRUG:
        y = keys.index(key)
        rows = drug[drug.label == DRUG[key]]
        assert len(rows)==23
        _, color, _ = STYLE[key]
        means = []
        for metric, dy, marker, filled in [('cov_p_075', .12, 'o', False), ('cov_r_075', -.12, 's', True)]:
            values = rows[metric].to_numpy() * 100
            assert np.isfinite(values).all()
            mean = values.mean(); lo, hi = np.quantile(values[rng.integers(23, size=(10000,23))].mean(1), [.025, .975])
            means.append(mean)
            a.errorbar(mean, y+dy, xerr=[[mean-lo],[hi-mean]], fmt=marker, color=color,
                       mfc=color if filled else 'white', ms=4.5, elinewidth=.8, capsize=2, mew=.9)
            exports.append(dict(key=key, metric=metric, mean=mean, ci_low=lo, ci_high=hi, bootstrap_n=10000, seed=20260924))
        a.plot(means, [y+.12,y-.12], lw=.6, color=color, alpha=.5)
    a.set_yticks(range(len(keys)), [STYLE[k][0] for k in keys])
    a.set(xlim=(15,103), ylim=(5.7,-.65), xlabel='Coverage (%)')
    style_axis(a, 'A', 'Observed-reference coverage', grid='x')
    fig.legend(handles=[Line2D([], [], marker='s', ls='none', color='#333333', label='Recall (COV-R)'),
                      Line2D([], [], marker='o', ls='none', color='#333333', mfc='white', label='Precision (COV-P)')],
             loc='lower center', bbox_to_anchor=(.58, .08), ncol=2, frameon=False, fontsize=7)
    q = drug[drug.label==QWEN].set_index('smiles')
    l = drug[drug.label==DRUG['loqi_raw']].set_index('smiles').reindex(q.index)
    x, y = 100*l.cov_r_075, 100*q.cov_r_075
    b.plot([0,100],[0,100], ls='--', lw=.8, color='#A8ADB5', zorder=1)
    points = pd.DataFrame({'x': x, 'y': y}).groupby(['x', 'y']).size().reset_index(name='n')
    assert points.n.sum() == 23
    b.scatter(points.x, points.y, s=22 + 12*(points.n-1), c='#A2A9B1',
              edgecolors='white', linewidth=.4, zorder=3)
    for point in points[points.n > 1].itertuples():
        b.annotate(f'{point.n} molecules', (point.x, point.y), xytext=(-62, -18),
                   textcoords='offset points', fontsize=6.5,
                   arrowprops={'arrowstyle': '-', 'lw': .6, 'color': '#6B7280'})
    for name, offset in [('Imatinib',(10,8)), ('Actinonin',(-43,13))]:
        select = q.name.str.lower().eq(name.lower()); assert select.sum()==1
        xx,yy=float(x[select].iloc[0]),float(y[select].iloc[0])
        b.scatter([xx],[yy], s=39, c=STYLE[QWEN if yy>xx else 'loqi_raw'][1], edgecolors='white', linewidth=.5, zorder=4)
        b.annotate(name,(xx,yy),xytext=offset,textcoords='offset points',fontsize=7,
                   arrowprops={'arrowstyle':'-', 'lw':.6,'color':'#6B7280'})
    b.set(xlim=(-4,104),ylim=(-4,104),xlabel='LoQI reference coverage (%)',ylabel='Qwen reference coverage (%)',aspect='equal')
    b.set_xticks([0,25,50,75,100]);b.set_yticks([0,25,50,75,100])
    style_axis(b,'B','Molecule-level comparison')
    fig.text(.5,.025,'Exploratory · 23 molecules · supplied pools without a common PoseBusters filter',ha='center',fontsize=7,color='#555555')
    pd.DataFrame(exports).to_csv(DATA/'drug_coverage_intervals.csv',index=False)
    save(fig,'multireference-coverage')


def historical_energy(energy):
    # This is deliberately NOT included in the manuscript: checkpoint identities
    # and corrected pool coverage differ from the main panel.
    fig,(a,b)=plt.subplots(1,2,figsize=(7,3.7))
    fig.subplots_adjust(left=.09,right=.98,top=.77,bottom=.25,wspace=.30)
    selection=[('qwen_1p7b_fsq_bigdata_pretrain_fixed','Legacy Qwen 1.7B FSQ',STYLE[QWEN][1],'o'),
               ('loqi_raw_fixed','Historical LoQI',STYLE['loqi_raw'][1],'s'),
               ('rdkit_random_raw_fixed','Historical RDKit raw',STYLE['rdkit_random_raw'][1],'X'),
               ('chembl3d_gt_pb','Historical ChEMBL3D-PB',STYLE['chembl3d_gt_pb'][1],'h')]
    used=[]
    for method,label,color,marker in selection:
        rows=energy[(energy.ligand_set=='core') & (energy.method==method)].set_index('energy_window')
        windows=['deltaE_5','deltaE_10','deltaE_20','deltaE_50','all']
        rows=rows.loc[windows]
        # Extend conditional published aggregate to the complete 94-entry universe.
        hit=100*rows.hit_0p75_window*rows.n_ligands/94
        count=rows.mean_n_window_confs*rows.n_ligands/94
        a.plot(range(5),hit,color=color,marker=marker,label=label)
        b.plot(range(5),count,color=color,marker=marker)
        for w,h,c in zip(windows,hit,count):
            used.append(dict(method=method,window=w,hit_full_core=h,mean_retained_full_core=c))
    for ax in (a,b):
        ax.set_xticks(range(5),['5','10','20','50','All'])
        ax.set_xlabel('Within-pool energy window (kcal/mol)')
    a.set(ylabel='Recovered targets at 0.75 Å (%)',ylim=(0,103))
    b.set(ylabel='Mean retained conformers',ylim=(0,1050))
    style_axis(a,'A','Historical recovery · core')
    style_axis(b,'B','Historical retained counts · core')
    fig.suptitle('Historical sidecar — not the corrected main comparison',fontsize=10,color='#8C3D26',y=.98)
    fig.legend(*a.get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,.015),fontsize=7)
    pd.DataFrame(used).to_csv(DATA/'historical_energy_plotted.csv',index=False)
    save(fig,'energy-window-historical')


def main():
    global GALLERY
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh-data', action='store_true')
    parser.add_argument('--historical-energy', action='store_true',
                        help='Also render the separate historical diagnostic.')
    args = parser.parse_args()
    if args.refresh_data:
        refresh_data()
    OUT.mkdir(exist_ok=True)
    with PdfPages(OUT / 'main-figures.pdf') as gallery:
        GALLERY = gallery
        recovery(pd.read_csv(DATA / 'recovery-table-tiers.csv'))
        diversity(pd.read_csv(DATA / 'core-comparison.csv'))
        size_figure(pd.read_csv(DATA / 'core-strata.csv'))
        drug_figure(pd.read_csv(DATA / 'drug_molecule_metrics.csv'))
    GALLERY = None
    if args.historical_energy:
        historical_energy(pd.read_csv(DATA / 'historical_energy_summary.csv'))
    print('Wrote four main figures (PDF, SVG, PNG) and the combined gallery.')


if __name__ == '__main__':
    main()
