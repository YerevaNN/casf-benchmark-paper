"""Render manuscript Figures 3, 4 and B1 from archived core/ref measurements."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parent
PAPER = OUT.parents[1]
DATA = PAPER / 'results-data/three-dataset-tables'
STYLE = json.loads((PAPER / 'sampling-data/palette.json').read_text())
KEYS = {value[0]: key for key, value in STYLE.items()}
TITLES = {'core': 'CASF core · 94 entries', 'ref': 'CASF ref · 1,236 entries'}
INK, MUTED = '#29323C', '#626C78'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10,
    'axes.titleweight': 'medium', 'axes.labelsize': 9, 'xtick.labelsize': 8,
    'ytick.labelsize': 8, 'axes.edgecolor': '#8B929A', 'axes.labelcolor': INK,
    'text.color': INK, 'xtick.color': '#535D68', 'ytick.color': '#535D68',
    'axes.linewidth': .65, 'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white',
    'figure.facecolor': 'white',
})
inputs = [DATA / 'table3.csv', PAPER / 'energy-data/summary.csv',
          PAPER / 'energy-data/protocol.json', PAPER / 'sampling-data/palette.json',
          PAPER / 'sampling-data/recovery_curves.csv', PAPER / 'sampling-data/recovery80_frontier.csv',
          PAPER / 'sampling-data/reference.json', PAPER / 'sampling-data/display.json']
table3 = pd.read_csv(inputs[0])
energy = pd.read_csv(inputs[1])
ref_results = OUT / "ref"
ref_energy = pd.read_csv(ref_results / "summary.csv")
inputs.extend(ref_results/name for name in ["summary.csv", "provenance.json", "recovery_curves.csv", "recovery80_frontier.csv", "full-pool-endpoints.csv", "reference.json"])
figures = []


def save(fig, name):
    for ext in ['pdf', 'svg', 'png']:
        path = PAPER / 'figures' / f'{name}.{ext}'
        fig.savefig(path, dpi=220)
        if ext == 'svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines()) + '\n')
    figures.append(fig)


def row_axis(ax, methods, show_names=True):
    ax.set_ylim(len(methods)-.4, -.6)
    ax.set_yticks(range(len(methods)), methods if show_names else [])
    ax.tick_params(axis='y', length=0)
    ax.axhline(len(methods)-1.5, color='#7B838D', lw=.7, ls=(0, (4, 3)))
    ax.grid(axis='x', color='#E7EBEE', linewidth=.6)
    ax.set_axisbelow(True)


# Figure 3: identical energy protocol and panel scales for both datasets.
fig, axes = plt.subplots(2, 2, figsize=(8.5, 9.4))
fig.subplots_adjust(left=.215, right=.965, top=.865, bottom=.16, hspace=.55, wspace=.28)
methods = energy.label.tolist()
y = np.arange(len(methods))
energy_points = []
for i, cohort in enumerate(['core', 'ref']):
    for j, metric in enumerate(['clusters', 'energy']):
        ax = axes[i, j]
        title = 'Geometric diversity' if j == 0 else 'Energy dispersion'
        ax.set_title(f'{chr(65+2*i+j)}   {TITLES[cohort]}\n{title}', loc='left', pad=12)
        row_axis(ax, methods, show_names=j == 0)
        ax.set_xlim(0, 66 if j == 0 else 49)
        ax.set_xlabel('Mean clusters at 0.5 Å' if j == 0 else 'Median energy SD (kcal/mol)', labelpad=7)
        for k, method in enumerate(methods):
            cohort_energy = energy if cohort == 'core' else ref_energy
            r = cohort_energy.loc[cohort_energy.label.eq(method)].iloc[0]
            value = float(r.mean_clusters_0p5 if j == 0 else r.median_energy_std)
            color = STYLE[KEYS[method]][1]
            ax.barh(k, value, color=color, height=.58, zorder=3)
            ax.text(value + (.9 if j == 0 else .65), k, f'{value:.1f}', va='center', fontsize=8)
            energy_points.append(dict(cohort=cohort, method=method, metric=metric, value=value))
fig.suptitle('ChEMBL-count target · geometry and energy spread', x=.215, ha='left', y=.974, fontsize=11)
fig.text(.5, .045, 'Both datasets: common MMFF94s protocol, rebuilt hydrogens, no coordinate optimization.\nEnergy dispersion is the median within-ensemble population SD across molecules.',
         ha='center', fontsize=8, color=MUTED, linespacing=1.5)
save(fig, 'figure-energy-dispersion')
pd.DataFrame(energy_points).to_csv(OUT/'figure3-measured-values.csv', index=False)

# Figure 4: every endpoint is available in current Table 3.
fixed = table3.loc[table3.cohort.eq('core') & table3.budget.eq('fixed')]
methods = fixed.sort_values('distance', ascending=False).method.tolist() + ['ChEMBL3D-PB']
fig, axes = plt.subplots(2, 2, figsize=(8.5, 9.4))
fig.subplots_adjust(left=.215, right=.965, top=.865, bottom=.18, hspace=.55, wspace=.28)
for i, cohort in enumerate(['core', 'ref']):
    for j, metric in enumerate(['recovery', 'clusters']):
        ax = axes[i, j]
        title = 'Experimental recovery' if j == 0 else 'Geometric diversity'
        ax.set_title(f'{chr(65+2*i+j)}   {TITLES[cohort]}\n{title}', loc='left', pad=12)
        row_axis(ax, methods, show_names=j == 0)
        ax.set_xlim((58, 100) if j == 0 else (0, 330))
        ax.set_xticks([60, 70, 80, 90, 100] if j == 0 else [0, 100, 200, 300])
        ax.set_xlabel('Recovery at 0.75 Å (%)' if j == 0 else 'Mean clusters at 0.5 Å', labelpad=7)
        for k, method in enumerate(methods):
            ax.axhspan(k-.42, k+.42, color='#F5F7F9' if k%2 == 0 else 'white', zorder=-1)
            r = table3.loc[table3.cohort.eq(cohort) & table3.method.eq(method)].set_index('budget')
            x0 = r.loc['chembl_count', metric]
            color = STYLE[KEYS[method]][1]
            if method == 'ChEMBL3D-PB':
                ax.scatter(x0, k, s=50, marker='h', color=color, edgecolor='white', linewidth=.65, zorder=5)
            else:
                x1 = r.loc['fixed', metric]
                ax.plot([x0, x1], [k, k], color=color, lw=1.8, alpha=.6, zorder=2)
                ax.scatter(x0, k, s=33, facecolor='white', edgecolor=color, linewidth=1.3, zorder=3)
                ax.scatter(x1, k, s=37, color=color, edgecolor='white', linewidth=.5, zorder=4)
handles = [Line2D([], [], marker='o', color='#697581', mfc='white', ls='none', label='ChEMBL-count target'),
           Line2D([], [], marker='o', color='#697581', ls='none', label='1,000-candidate target'),
           Line2D([], [], marker='h', color='#333A43', ls='none', label='Stored ChEMBL3D-PB')]
fig.suptitle('Candidate targets · recovery and geometric diversity', x=.215, ha='left', y=.974, fontsize=11)
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.55, .07), ncol=3,
           frameon=False, fontsize=8, columnspacing=1.1, handletextpad=.5)
fig.text(.5, .032, 'Lines join two evaluated targets; retained counts and computing costs are not matched.',
         ha='center', fontsize=7.7, color=MUTED)
save(fig, 'figure-3-sampling-budget')
table3.to_csv(OUT/'figure4-source-values.csv', index=False)

# Figure B1: core curves are archived; ref full-pool endpoints are independently available.
curves = pd.read_csv(PAPER/'sampling-data/recovery_curves.csv')
frontier = pd.read_csv(PAPER/'sampling-data/recovery80_frontier.csv')
core_ref = json.loads((PAPER/'sampling-data/reference.json').read_text())
display = json.loads((PAPER/'sampling-data/display.json').read_text())
methods = display['methods']
# Independently checked against the published per-entry database in the calculation audit.
ref_ends = pd.read_csv(ref_results/'full-pool-endpoints.csv').set_index('method')
for method, row in ref_ends.iterrows():
    match = table3.loc[table3.cohort.eq('ref') &
                       table3.method.eq(STYLE[method.removesuffix('_fixed')][0]) &
                       table3.budget.eq('chembl_count' if method == 'chembl3d_gt_pb' else 'fixed')]
    assert len(match) == 1 and np.isclose(row.recovery, match.recovery.iloc[0])
fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.2))
fig.subplots_adjust(left=.075, right=.85, top=.86, bottom=.17, hspace=.61, wspace=.87)
short = {methods[0]: 'Qwen', methods[1]: 'LoQI', methods[2]: 'FlowR', methods[3]: 'NExT-Mol'}

def annotate_endpoints(ax, values, precision, unit, min_gap):
    # Place labels in vertical order, preserving genuine endpoint locations.
    ordered = sorted(values.items(), key=lambda pair: pair[1])
    label_positions = []
    for method, value in ordered:
        label_positions.append(max(value, label_positions[-1]+min_gap) if label_positions else value)
    shift = max(0, label_positions[-1] - (ax.get_ylim()[1] - min_gap/2))
    for (method, value), position in zip(ordered, label_positions):
        color = STYLE[method.removesuffix('_fixed')][1]
        ax.plot(1000, value, marker='o', ms=3.5, color=color, clip_on=False)
        ax.annotate(f'{short[method]}  {value:.{precision}f}{unit}', xy=(1000,value),
                    xytext=(1.045,position-shift), textcoords=ax.get_yaxis_transform(),
                    ha='left', va='center', fontsize=8, color=color, annotation_clip=False,
                    arrowprops=dict(arrowstyle='-',color=color,lw=.7))

for i, cohort in enumerate(['core', 'ref']):
    a, b = axes[i]
    for ax in (a,b):
        ax.set_xlim(10,1000)
        ax.set_xticks([10,100,250,500,1000], ['10','100','250','500','1,000'])
        ax.grid(color='#E7EBEE',lw=.6)
        ax.set_axisbelow(True)
        ax.set_xlabel('Candidate budget')
    a.set_ylim(50,95); a.set_yticks([50,60,70,80,90,95]); a.set_ylabel('Expected recovery (%)')
    b.set_ylim(0,1.5); b.set_yticks([0,.3,.6,.9,1.2,1.5]); b.set_ylabel('RMSD threshold (Å)')
    a.set_title(f'{chr(65+2*i)}   {TITLES[cohort]}\nRecovery at 0.75 Å', loc='left',pad=12)
    b.set_title(f'{chr(66+2*i)}   {TITLES[cohort]}\nThreshold for 80% recovery',loc='left',pad=12)
    ends_a, ends_b = {}, {}
    for method in methods:
        color = STYLE[method.removesuffix('_fixed')][1]
        cohort_curves = curves if cohort == 'core' else pd.read_csv(ref_results/'recovery_curves.csv')
        cohort_frontier = frontier if cohort == 'core' else pd.read_csv(ref_results/'recovery80_frontier.csv')
        c = cohort_curves.loc[cohort_curves.method.eq(method) & cohort_curves.threshold.eq(.75)]
        f = cohort_frontier.loc[cohort_frontier.method.eq(method)]
        a.plot(c.budget,100*c.recovery,color=color,lw=1.9)
        b.plot(f.budget,f.rmsd_threshold,color=color,lw=1.9)
        ends_a[method] = 100*c.loc[c.budget.eq(1000), 'recovery'].iloc[0]
        ends_b[method] = f.loc[f.budget.eq(1000), 'rmsd_threshold'].iloc[0]
        if cohort == 'ref':
            assert np.isclose(ends_a[method], ref_ends.loc[method, 'recovery'])
            assert np.isclose(ends_b[method], ref_ends.loc[method, 'rmsd_threshold'])
    reference = [100*core_ref['hit_rates']['0.75'],core_ref['rmsd_threshold']] if cohort=='core' else ref_ends.loc['chembl3d_gt_pb',['recovery','rmsd_threshold']].tolist()
    for ax, value in zip((a,b),reference):
        ax.axhline(value,color=STYLE['chembl3d_gt_pb'][1],ls=(0,(4,3)),lw=1.2)
    annotate_endpoints(a,ends_a,1,'%',3.4)
    annotate_endpoints(b,ends_b,3,' Å',.12)
handles = [Line2D([],[],color=STYLE[m.removesuffix('_fixed')][1],lw=1.9,label=STYLE[m.removesuffix('_fixed')][0]) for m in methods]
handles.append(Line2D([],[],color=STYLE['chembl3d_gt_pb'][1],ls=(0,(4,3)),lw=1.2,label='Stored ChEMBL3D-PB'))
fig.suptitle('Sampling budget, recovery, and matching tolerance',x=.075,ha='left',y=.98,fontsize=11.5)
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.065),ncol=3,frameon=False,fontsize=8,columnspacing=1.5)
fig.text(.5,.025,'Rejected candidates consume budget. Dashed lines show each dataset’s full stored ChEMBL3D-PB ensemble.',ha='center',fontsize=7.5,color=MUTED)
save(fig,'figure-budget-combined')
for fig in figures:
    plt.close(fig)
(OUT/'provenance.json').write_text(json.dumps({
    'analysis_commit':'5090d110671731fefbe31c81fb65258c178e6304',
    'release':'https://github.com/YerevaNN/casf-benchmark/releases/tag/dashboard-data-ref-protocol-v1',
    'source_hashes':{str(p.relative_to(PAPER)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    'figure3':'Core and ref common-hydrogen MMFF94s energies; matched geometric summaries.',
    'figure4':'All core/ref endpoints from current Table 3; same method order and scales.',
    'figureB1':'Both datasets: exact pre-PoseBusters finite-pool curves; rejected candidates consume budget. Ref endpoint and threshold checks agree with archived per-entry minima.',
    'qwen':'Original 1.7B FSQ step47023 in all panels.',
    'missing_inputs':[]
},indent=2)+'\n')
print('Completed Figures 3, 4 and B1 with calculated core/ref energies and exact sampling curves.')
