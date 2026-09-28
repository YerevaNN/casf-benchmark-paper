"""Recreate manuscript figures from the public dashboard release.

Run from the project environment: python manuscript/build_figures.py
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
        'method_selection': 'Main curve panel excludes minimized classical variants and additional Qwen checkpoints; full main comparison remains in Table 1.',
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
    fig.savefig(OUT / f'{name}.svg')
    fig.savefig(OUT / f'{name}.png', dpi=450)
    plt.close(fig)


def get_row(summary, cohort, key, tier='fixed'):
    if key == 'chembl3d_gt_pb':
        tier = 'reference'
    selected = summary[(summary.cohort == cohort) & (summary.key == key) & (summary.tier == tier)]
    assert len(selected) == 1
    return selected.iloc[0]


def recovery(summary):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7, 4.35), gridspec_kw={'width_ratios': [1.05, 1.0]})
    fig.subplots_adjust(left=.085, right=.98, top=.86, bottom=.25, wspace=.60)
    for key in reversed(KEYS):
        _, color, marker = STYLE[key]
        row = get_row(summary, 'core', key)
        a.plot([v for _, v in THRESHOLDS], [row['hit_' + t] for t, _ in THRESHOLDS],
               marker=marker, color=color, ls='--' if key == 'chembl3d_gt_pb' else '-',
               lw=1.8 if key == QWEN else 1.1, zorder=4 if key == QWEN else 2)
    a.axvline(.75, color='#B6BDC5', ls=':', lw=.8)
    a.set(xlabel='Matching RMSD cutoff (Å)', ylabel='Recovered targets (%)', ylim=(15, 103), xlim=(.18, 2.06))
    a.set_xticks([.25, .5, .75, 1.0, 1.5, 2.0], ['0.25', '0.5', '0.75', '1.0', '1.5', '2.0'])
    style_axis(a, 'A', 'Recovery thresholds · core (n = 94)')
    budget = KEYS[:-1]
    for y, key in enumerate(budget):
        name, color, _ = STYLE[key]
        lo = get_row(summary, 'core', key, 'chembl_count').hit_0p75
        hi = get_row(summary, 'core', key).hit_0p75
        b.plot([lo, hi], [y, y], color=color, lw=1.6)
        b.plot(lo, y, 'o', color=color, mfc='white', ms=5, mew=1.1)
        b.plot(hi, y, 's', color=color, ms=4.5)
    ref = get_row(summary, 'core', 'chembl3d_gt_pb').hit_0p75
    b.axvline(ref, color=STYLE['chembl3d_gt_pb'][1], ls='--', lw=1)
    b.set_yticks(range(len(budget)), [STYLE[k][0] for k in budget])
    b.set_ylim(len(budget)-.4, -.7)
    b.set(xlabel='Recovered targets at 0.75 Å (%)', xlim=(50, 100))
    style_axis(b, 'B', 'Candidate-target comparison', grid='x')
    b.legend(handles=[Line2D([], [], marker='o', color='#333333', mfc='white', ls='none', label='ChEMBL-count target'),
                      Line2D([], [], marker='s', color='#333333', ls='none', label='1,000-candidate target')],
             loc='upper left', bbox_to_anchor=(-.04, 1.02), frameon=False, fontsize=6.4, borderpad=0)
    # Extra headroom keeps the marker key separate from the first method row.
    b.set_ylim(len(budget)-.4, -1.8)
    legend(fig)
    save(fig, 'recovery-thresholds-budget')


def diversity(summary):
    radii = [(tag, radius) for tag, radius in RADII if radius <= 2.0]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7, 4.1))
    fig.subplots_adjust(left=.09, right=.975, top=.85, bottom=.26, wspace=.30)
    for key in reversed(KEYS):
        name, color, marker = STYLE[key]
        row = get_row(summary, 'core', key)
        a.plot([v for _, v in radii], [row['clusters_' + t] for t, _ in radii], marker=marker,
               color=color, ls='--' if key == 'chembl3d_gt_pb' else '-', lw=1.8 if key == QWEN else 1.1)
        b.scatter(row.clusters_1p0, row.hit_0p75, marker=marker, s=36, c=color, edgecolors='white', linewidths=.4, zorder=3)
        offsets = {QWEN:(5, 7), 'loqi_raw':(-4, 8), 'flowr_raw':(-31, 9), 'torsion_raw':(5, -12), 'chembl3d_gt_pb':(4,-14)}
        if key in offsets:
            b.annotate(name, (row.clusters_1p0, row.hit_0p75), xytext=offsets[key], textcoords='offset points', fontsize=7, color=color)
    a.axvline(1, color='#B6BDC5', ls=':', lw=.8)
    a.set(xlabel='Clustering RMSD radius (Å)', ylabel='Mean cluster count', xlim=(.4, 2.1), ylim=(0, None))
    a.set_xticks([radius for _, radius in radii])
    b.set(xlabel='Mean cluster count at 1.0 Å', ylabel='Recovered targets at 0.75 Å (%)', xlim=(0, 125), ylim=(72, 99))
    style_axis(a, 'A', 'Geometric resolution · core')
    style_axis(b, 'B', 'Diversity and recovery · core')
    legend(fig)
    save(fig, 'diversity-recovery')


def size_figure(strata):
    fig, axes = plt.subplots(2, 2, figsize=(7, 5.7), sharey=True)
    fig.subplots_adjust(left=.09, right=.98, top=.92, bottom=.17, wspace=.20, hspace=.43)
    for i, cohort in enumerate(['core', 'ref']):
        for j, descriptor in enumerate(['size', 'rotors']):
            ax = axes[i, j]
            example = strata[(strata.cohort == cohort) & (strata.key == QWEN) & (strata.descriptor == descriptor)]
            labels = example.group.tolist()
            for key in reversed(FOCUS):
                rows = strata[(strata.cohort == cohort) & (strata.key == key) & (strata.descriptor == descriptor)].set_index('group').loc[labels]
                _, color, marker = STYLE[key]
                ax.plot(range(4), rows.hit, marker=marker, color=color, ls='--' if key=='chembl3d_gt_pb' else '-', lw=1.6)
            ax.set_xticks(range(4), [f'{g}\nn = {n}' for g,n in zip(labels,example.n)])
            ax.set(ylim=(-3, 106), xlim=(-.2, 3.2), xlabel='Heavy atoms' if descriptor=='size' else 'Rotatable bonds')
            if j == 0:
                ax.set_ylabel('Recovered targets at 0.75 Å (%)')
            title = f'{cohort.capitalize()} · '+('molecular size' if descriptor=='size' else 'flexibility')
            style_axis(ax, 'ABCD'[i*2+j], title)
            if i == 0:
                ax.axvspan(2.65, 3.2, color='#F1F3F5', zorder=0)
    legend(fig, FOCUS, ncol=4, y=.025)
    save(fig, 'size-flexibility')


def drug_figure(drug):
    rng = np.random.default_rng(20260924)
    keys = list(DRUG)
    fig, (a, b) = plt.subplots(1, 2, figsize=(7, 4.05), gridspec_kw={'width_ratios':[1,1]})
    fig.subplots_adjust(left=.19, right=.98, top=.85, bottom=.24, wspace=.37)
    exports = []
    for y,key in enumerate(keys):
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
    b.scatter(x,y, s=22, c='#A2A9B1', edgecolors='white', linewidth=.4, zorder=3)
    for name, offset in [('Imatinib',(10,8)), ('Actinonin',(-43,13))]:
        select = q.name.str.lower().eq(name.lower()); assert select.sum()==1
        xx,yy=float(x[select].iloc[0]),float(y[select].iloc[0])
        b.scatter([xx],[yy], s=39, c=STYLE[QWEN if yy>xx else 'loqi_raw'][1], edgecolors='white', linewidth=.5, zorder=4)
        b.annotate(name,(xx,yy),xytext=offset,textcoords='offset points',fontsize=7,
                   arrowprops={'arrowstyle':'-', 'lw':.6,'color':'#6B7280'})
    b.set(xlim=(-4,104),ylim=(-4,104),xlabel='LoQI reference coverage (%)',ylabel='Qwen reference coverage (%)',aspect='equal')
    b.set_xticks([0,25,50,75,100]);b.set_yticks([0,25,50,75,100])
    style_axis(b,'B','Molecule-level comparison')
    fig.text(.5,.025,'23 molecules · supplied pools without a common PoseBusters filter',ha='center',fontsize=7,color='#555555')
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
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh-data',action='store_true')
    args=parser.parse_args()
    if args.refresh_data:
        refresh_data()
    summary=pd.read_csv(DATA/'casf_summary.csv')
    recovery(summary)
    diversity(summary)
    size_figure(pd.read_csv(DATA/'casf_strata.csv'))
    drug_figure(pd.read_csv(DATA/'drug_molecule_metrics.csv'))
    historical_energy(pd.read_csv(DATA/'historical_energy_summary.csv'))
    print('Wrote four manuscript figures and one separate historical diagnostic (PDF, SVG, PNG).')


if __name__=='__main__':
    main()
