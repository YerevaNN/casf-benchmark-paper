"""Render the optional figure-review appendix from archived corrected data.

Run: python jcim-overleaf/build_candidate_figures.py
No database access, manuscript compilation, or changes to the four main figures.
"""
from pathlib import Path
import json
import hashlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from build_figures import ROOT, DATA, PANEL, STYLE, QWEN, DRUG

OUT = ROOT / 'figures/candidates'
REF = 'chembl3d_gt_pb'
KEYS = [REF] + PANEL
DK = [k for k in PANEL if k in DRUG]
SEED = 20260928
NBOOT = 10000
CATALOG = []
GALLERY = None


def label(key):
    return STYLE[key][0]


def finish(fig, slug, title, caption, source, drug=False):
    number = len(CATALOG) + 1
    code = f'C{number:02d}'
    fig.suptitle(f'{code}  |  {title}', fontsize=11, x=.02, ha='left')
    status = ('Exploratory drug set · 23 molecules · supplied pools without common PB filtering'
              if drug else 'CASF core · 94 ligands · recovery cutoff 0.75 Å')
    fig.text(.02, .012, status, fontsize=7, color='#555555')
    fig.tight_layout(rect=(0, .045, 1, .94))
    stem = f'{code.lower()}-{slug}'
    fig.savefig(OUT / f'{stem}.pdf')
    fig.savefig(OUT / f'{stem}.png', dpi=220)
    fig.savefig(OUT / f'{stem}.svg')
    svg = OUT / f'{stem}.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
    GALLERY.savefig(fig)
    plt.close(fig)
    CATALOG.append(dict(id=code, stem=stem, title=title, caption=caption, source=source))


def method_axis(ax, keys=KEYS):
    ax.set_yticks(range(len(keys)), [label(k) for k in keys])
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=.15)
    ax.set_axisbelow(True)


def annotated_heatmap(ax, values, rows, cols, vmax=100, cmap='Blues', fmt='.0f'):
    im = ax.imshow(values, aspect='auto', vmin=0, vmax=vmax, cmap=cmap)
    ax.set_yticks(range(len(rows)), rows)
    ax.set_xticks(range(len(cols)), cols, rotation=45, ha='right')
    for i, j in np.ndindex(values.shape):
        value = values[i, j]
        ax.text(j, i, format(value, fmt), ha='center', va='center', fontsize=6.5,
                color='white' if value > .58 * vmax else '#222222')
    return im


def bootstrap(delta, indices):
    means = delta[indices].mean(axis=1)
    return delta.mean(), *np.quantile(means, [.025, .975])


def core_figures(entries, tiers):
    fixed = entries[entries.method.str.endswith('_fixed') | entries.method.eq(REF)].copy()
    fixed['key'] = fixed.method.str.replace(r'_fixed$', '', regex=True)
    hits = fixed.pivot(index='mol_id', columns='key', values='casf_hit_0p75').reindex(columns=KEYS).fillna(0)
    assert hits.shape == (94, 11) and hits.isin([0, 1]).all().all()
    expected = tiers[(tiers.cohort == 'core') & tiers.tier.isin(['fixed', 'reference'])].set_index('key')
    np.testing.assert_allclose(hits.mean() * 100, expected.reindex(KEYS).hit_0p75)
    base = fixed[fixed.key.eq(REF)].set_index('mol_id').reindex(hits.index)
    mean = fixed.groupby('key').mean(numeric_only=True).reindex(KEYS)
    rng = np.random.default_rng(SEED)
    indices = rng.integers(0, len(hits), (NBOOT, len(hits)))
    intervals = []
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for y, key in enumerate(PANEL):
        delta, low, high = bootstrap((hits[key] - hits[REF]).to_numpy() * 100, indices)
        ax.errorbar(delta, y, xerr=[[delta-low], [high-delta]], fmt=STYLE[key][2], color=STYLE[key][1], capsize=3)
        intervals.append(dict(key=key, delta=delta, ci_low=low, ci_high=high))
    method_axis(ax, PANEL)
    ax.axvline(0, color='#444444', ls='--', lw=.8)
    ax.set_xlabel('Recovery difference from ChEMBL3D-PB (percentage points)')
    pd.DataFrame(intervals).to_csv(OUT / 'core-paired-intervals.csv', index=False)
    finish(fig, 'paired-recovery', 'How much recovery improves over the stored pool',
           r'Paired differences in core recovery between each 1,000-candidate pipeline and the finite stored ChEMBL3D-PB pool. Error bars are exploratory 95\% percentile intervals from 10,000 paired bootstrap resamples of the 94 core entries (94 distinct parent compounds; seed 20260928). Missing outputs are misses. The comparison does not match retained count or computing cost, and intervals are not adjusted for multiple comparisons or checkpoint selection.', 'HTML §§3–4; core-entries.csv')

    fig, ax = plt.subplots(figsize=(8, 4.5))
    max_lost, max_new = 0, 0
    for y, key in enumerate(PANEL):
        new = int(((hits[key] == 1) & (hits[REF] == 0)).sum())
        lost = int(((hits[key] == 0) & (hits[REF] == 1)).sum())
        max_lost, max_new = max(max_lost, lost), max(max_new, new)
        ax.barh(y, new, color=STYLE[key][1], height=.65)
        ax.barh(y, -lost, color=STYLE[key][1], alpha=.4, height=.65)
        ax.text(new+.25, y, str(new), va='center', fontsize=7)
        ax.text(-lost-.25, y, str(lost), va='center', ha='right', fontsize=7)
    method_axis(ax, PANEL); ax.axvline(0, color='black', lw=.7)
    ax.set_xlim(-max_lost-2, max_new+2); ax.set_xlabel('Lost recoveries  ←  Number of ligands  →  Additional recoveries')
    finish(fig, 'gained-lost', 'Additional recoveries and lost recoveries',
           r'Numbers of core ligands recovered only by a generated pool (right) or only by ChEMBL3D-PB (left). Shared hits and shared misses are omitted. Generated pools use the 1,000-candidate target and PoseBusters filtering; the reference is the finite stored pool. A positive net difference need not imply improvement on every ligand.', 'HTML §3; core-entries.csv')

    order = hits.assign(total=hits.sum(axis=1)).sort_values(['total'] + KEYS, kind='stable').index
    fig, ax = plt.subplots(figsize=(11, 4.7))
    ax.imshow(hits.loc[order].T, aspect='auto', cmap=ListedColormap(['#ececec', '#2171b5']), vmin=0, vmax=1)
    short_ids = [value.split('_')[0] for value in order]
    assert len(set(short_ids)) == 94
    ax.set_yticks(range(11), [label(k) for k in KEYS]); ax.set_xticks(range(94), short_ids, rotation=90, fontsize=5)
    ax.set_xlabel('Core ligand IDs, ordered from fewest to most methods recovering them')
    ax.legend(handles=[Line2D([], [], marker='s', ls='', color='#2171b5', label='Recovered'), Line2D([], [], marker='s', ls='', color='#dddddd', label='Missed')], loc='upper left', bbox_to_anchor=(0, 1.13), ncol=2)
    hits.loc[order].to_csv(OUT / 'core-recovery-matrix.csv')
    finish(fig, 'ligand-recovery', 'Which ligands are recovered by which methods?',
           r'Binary recovery for every core ligand and each of the ten generation pipelines plus ChEMBL3D-PB. Blue cells indicate a PoseBusters-passing conformer within 0.75~\AA{}; gray cells indicate a miss, including missing output. Columns are sorted by the number of methods recovering each ligand. Generated pools use the 1,000-candidate target. This is a descriptive map of correlated outcomes, not an independent set of method trials.', 'HTML §§3,6; core-entries.csv')

    values = np.array([[((hits[a] == 1) & (hits[b] == 0)).sum() for b in KEYS] for a in KEYS])
    fig, ax = plt.subplots(figsize=(8, 6.5))
    annotated_heatmap(ax, values, [label(k) for k in KEYS], [label(k) for k in KEYS], vmax=max(values.max(), 1))
    ax.set_xlabel('Method missing the ligand'); ax.set_ylabel('Method recovering the ligand')
    finish(fig, 'complementarity', 'Pairwise complementarity between generators',
           r'Each cell counts ligands recovered by the row method and missed by the column method. The matrix is directional; its diagonal is zero. Generated pools use 1,000 candidates before filtering, while ChEMBL3D-PB uses its stored pool. Counts identify complementary recovery patterns; they do not measure the performance or cost of a newly generated mixed ensemble.', 'HTML §§3,12; core-entries.csv')

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.7), sharey=True)
    for ax, tier, title in zip(axes, ['chembl_count', 'fixed'], ['ChEMBL-count target', '1,000-candidate target']):
        sel = tiers[(tiers.cohort == 'core') & (tiers.tier.eq(tier) | tiers.key.eq(REF))].set_index('key').reindex(KEYS)
        for y, key in enumerate(KEYS):
            a, b = sel.loc[key, ['mean_minimum_rmsd', 'mean_ensemble_median_rmsd']]
            ax.plot([a, b], [y, y], color=STYLE[key][1]); ax.scatter(a, y, c=STYLE[key][1], marker='s', s=24); ax.scatter(b, y, edgecolors=STYLE[key][1], facecolors='white', s=24)
        ax.set_title(title); ax.set_xlabel('RMSD to bound geometry (Å)'); ax.grid(axis='x', alpha=.15)
    axes[0].set_yticks(range(11), [label(k) for k in KEYS]); axes[0].invert_yaxis()
    axes[1].legend(handles=[Line2D([], [], marker='s', ls='', color='#444444', label='Mean minimum'), Line2D([], [], marker='o', ls='', color='#444444', markerfacecolor='white', label='Mean ensemble median')], loc='upper center', bbox_to_anchor=(.5, -.12), ncol=2)
    finish(fig, 'best-typical', 'Best available conformers versus typical conformers',
           r'Mean per-ligand minimum RMSD (filled squares) and mean per-ligand ensemble-median RMSD (open circles), at both candidate targets. These continuous distances are conditional on measurable ensembles, unlike recovery fractions, whose denominator includes all 94 core entries. ChEMBL3D-PB is the same stored pool in both panels. Connecting segments indicate the difference between two summaries, not confidence intervals or paired transformations.', 'HTML §3; recovery-table-tiers.csv')

    rmsd = fixed.pivot(index='mol_id', columns='key', values='casf_best_rmsd')
    pairs = [(REF, 'loqi_raw'), (REF, QWEN), ('loqi_raw', QWEN), ('flowr_raw', QWEN)]
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    for ax, (a, b) in zip(axes.flat, pairs):
        pts = rmsd[[a, b]].dropna(); lim = max(1, pts.max().max()) * 1.06
        ax.scatter(pts[a], pts[b], c=STYLE[b][1], s=13, alpha=.65)
        ax.plot([0, lim], [0, lim], color='#888888', lw=.8)
        ax.axvline(.75, color='#888888', ls=':', lw=.8); ax.axhline(.75, color='#888888', ls=':', lw=.8)
        ax.set(xlim=(0, lim), ylim=(0, lim), xlabel=label(a)+' minimum RMSD (Å)', ylabel=label(b)+' minimum RMSD (Å)')
        ax.text(.04, .94, f'{len(pts)} paired ligands', transform=ax.transAxes, va='top', fontsize=7)
    finish(fig, 'paired-rmsd', 'Per-ligand geometric agreement in selected comparisons',
           r'Minimum RMSD for paired measurable core ensembles in four selected comparisons. Each point is one ligand; the diagonal indicates equal distance and dotted lines mark the sole recovery cutoff, 0.75~\AA{}. Lower values are better. Points below the diagonal favor the vertical-axis method. Generated pools use 1,000 candidates before filtering. Panel-specific sample counts exclude unmeasurable pairs; they are not the recovery denominator.', 'HTML §§3,5; core-entries.csv')

    size = pd.cut(base.heavy_atoms, [0, 20, 30, 40, np.inf], right=False, labels=['<20', '20–29', '30–39', '≥40'])
    rotor = pd.cut(base.rotatable_bonds, [-1, 3, 6, 8, np.inf], labels=['0–3', '4–6', '7–8', '≥9'])
    counts = pd.crosstab(rotor, size, dropna=False)
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    annotated_heatmap(ax, counts.to_numpy(), counts.index.tolist(), counts.columns.tolist(), vmax=counts.to_numpy().max())
    ax.set_xlabel('Heavy atoms'); ax.set_ylabel('Rotatable bonds')
    finish(fig, 'cohort-composition', 'Where the core cohort has statistical support',
           r'Counts of the 94 core ligands cross-tabulated by heavy-atom and rotatable-bond groups using the same reference-derived descriptors for every method. Sparse large or flexible groups limit the stability of apparent method rankings. Size and flexibility are associated descriptors; this table does not isolate their independent effects.', 'HTML §§1,6; core-entries.csv')

    fig, axes = plt.subplots(1, 3, figsize=(10, 4.7), sharey=True)
    for ax, metric, title in zip(axes, ['greedy_clusters_1p0', 'clusters_per_100_1p0', 'cluster_entropy_1p0'], ['Mean clusters', 'Mean clusters per 100 conformers', 'Mean normalized cluster entropy']):
        for y, key in enumerate(KEYS): ax.scatter(mean.loc[key, metric], y, color=STYLE[key][1], marker=STYLE[key][2], s=35)
        ax.set_title(title, fontsize=8); ax.grid(axis='x', alpha=.15)
    axes[0].set_yticks(range(11), [label(k) for k in KEYS]); axes[0].invert_yaxis()
    finish(fig, 'diversity-profiles', 'Different summaries of geometric diversity',
           r'Per-ligand diversity summaries averaged over measurable core ensembles: cluster count, clusters per 100 conformers, and Shannon entropy divided by the logarithm of the number of occupied clusters (zero for a single cluster), at a 1.0~\AA{} clustering radius. This radius is distinct from the 0.75~\AA{} recovery cutoff. Generated pools use the 1,000-candidate target; ChEMBL3D-PB uses its finite pool. Normalizing cluster count does not fully remove sampling-size effects, and none of these metrics measures recovery by itself.', 'HTML §5; core-entries.csv')

    fig, ax = plt.subplots(figsize=(8, 4.8))
    for key in PANEL:
        ax.scatter(mean.loc[key, 'pb_pass_rate']*100, hits[key].mean()*100, color=STYLE[key][1], marker=STYLE[key][2], s=65, label=label(key))
    ax.set(xlabel='Mean per-ligand PoseBusters pass rate (%)', ylabel='Recovery at 0.75 Å (%)')
    ax.grid(alpha=.15); ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    finish(fig, 'validity-recovery', 'Physical validity and bound-conformer recovery',
           r'Mean per-ligand PoseBusters pass rate versus recovery for the ten generated core ensembles at the 1,000-candidate target. PoseBusters yield describes the analyzed input pools; it does not include every possible decoding or preprocessing rejection. Recovery requires a passing conformer and includes all 94 entries. High validity is a separate requirement from proximity to the observed bound conformation.', 'HTML §§3,9; core-entries.csv')

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    for ax, metric, title in zip(axes, ['casf_hit_0p75', 'greedy_clusters_1p0'], ['Recovery at 0.75 Å (%)', 'Mean clusters at 1.0 Å']):
        for offset, raw, minimized in [(-.06, 'rdkit_random_raw', 'rdkit_random_minimized'), (.06, 'torsion_raw', 'torsion_minimized')]:
            vals = mean.loc[[raw, minimized], metric].to_numpy() * (100 if metric == 'casf_hit_0p75' else 1)
            ax.plot(np.array([0, 1])+offset, vals, marker=STYLE[raw][2], color=STYLE[raw][1], label=label(raw).replace(' raw', ''))
            for x, value in zip(np.array([0, 1])+offset, vals): ax.annotate(f'{value:.1f}', (x, value), xytext=(5, 5), textcoords='offset points', fontsize=7)
        ax.set_xticks([0, 1], ['Raw pool', 'Minimized pool']); ax.set_xlim(-.3, 1.35); ax.set_ylabel(title); ax.grid(axis='y', alpha=.15)
    axes[0].legend()
    finish(fig, 'classical-pools', 'Classical raw and minimized pool comparisons',
           r'Recovery and geometric diversity of RDKit and torsion-based raw and minimized core pools at the 1,000-candidate target. These are separately generated endpoint pools, not paired before-and-after measurements on identical conformers. Lines identify related pipelines and do not establish a causal effect of minimization. Diversity is averaged over measurable ensembles; recovery includes all core ligands.', 'HTML §9; core-entries.csv')


def drug_figures(drug):
    drug = drug.copy(); drug['name'] = drug.name.str.strip()
    matrices = {metric: drug.pivot(index='smiles', columns='label', values=metric).reindex(columns=[DRUG[k] for k in DK]) for metric in ['cov_r_075', 'cov_p_075', 'pb_pass_rate', 'mat_r', 'mat_p']}
    assert all(m.shape == (23, 6) and m.notna().all().all() for m in matrices.values())
    names = drug.drop_duplicates('smiles').set_index('smiles').name
    order = names.sort_values().index
    def mat(metric): return matrices[metric].loc[order].to_numpy()
    labels = [label(k) for k in DK]
    common = r'All drug metrics describe supplied pools without a common PoseBusters mask; failed-RMSD handling also differs between imported evaluations. These are exploratory comparisons.'

    fig, axes = plt.subplots(1, 2, figsize=(10, 7.5))
    for ax, metric, title in zip(axes, ['cov_r_075', 'cov_p_075'], ['Reference coverage / recall (COV-R, %)', 'Observed-reference precision (COV-P, %)']):
        annotated_heatmap(ax, mat(metric)*100, names.loc[order].tolist(), labels)
        ax.set_title(title, fontsize=9)
    finish(fig, 'drug-coverage-matrix', 'Coverage and precision for every drug molecule',
           r'Molecule-level COV-R and COV-P at the strict 0.75~\AA{} cutoff for all six selected generators. The two panels share a 0--100\% color scale. COV-R measures the fraction of observed references reached; COV-P measures the fraction of generated conformers near an observed reference. Unmatched conformers are not necessarily biologically impossible. '+common, 'HTML §7; drug_molecule_metrics.csv', True)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    offsets = [(6, 10), (6, 4), (6, 12), (6, -13), (-8, -10), (6, -13)]
    for key, offset in zip(DK, offsets):
        row = drug[drug.label.eq(DRUG[key])]
        x, y = row.cov_r_075.mean()*100, row.cov_p_075.mean()*100
        ax.scatter(x, y, c=STYLE[key][1], marker=STYLE[key][2], s=55)
        ax.annotate(label(key), (x, y), xytext=offset, textcoords='offset points', fontsize=8, ha='right' if key == 'flowr_raw' else 'left')
    ax.set(xlabel='Mean reference coverage / recall (%)', ylabel='Mean observed-reference precision (%)', xlim=(75, 98), ylim=(30, 65)); ax.grid(alpha=.15)
    finish(fig, 'drug-tradeoff', 'Coverage and concentration define different preferences',
           r'Molecule-mean recall and precision at 0.75~\AA{} for the six selected generators. Each of the 23 molecules receives equal weight regardless of its reference count. Higher values indicate broader observed-reference coverage or greater concentration near the available observations; they do not measure binding affinity or conformer populations. '+common, 'HTML §§7,12; drug_molecule_metrics.csv', True)

    delta_r = (matrices['cov_r_075'][DRUG[QWEN]] - matrices['cov_r_075'][DRUG['loqi_raw']])*100
    ordered = delta_r.sort_values().index
    fig, axes = plt.subplots(1, 2, figsize=(9, 7), sharey=True)
    for ax, metric, title in zip(axes, ['cov_r_075', 'cov_p_075'], ['Recall difference (percentage points)', 'Precision difference (percentage points)']):
        delta = (matrices[metric][DRUG[QWEN]] - matrices[metric][DRUG['loqi_raw']]).loc[ordered]*100
        ax.hlines(range(23), 0, delta, color='#bbbbbb'); ax.scatter(delta, range(23), c=np.where(delta>=0, STYLE[QWEN][1], STYLE['loqi_raw'][1]), s=25)
        ax.axvline(0, color='#555555', lw=.8); ax.set_xlabel(title); ax.grid(axis='x', alpha=.15)
    axes[0].set_yticks(range(23), names.loc[ordered]); axes[0].invert_yaxis()
    finish(fig, 'drug-paired-molecules', 'Similar mean recall can hide opposing molecule-level outcomes',
           r'Paired molecule-level differences, Qwen 1.7B FSQ minus LoQI, in recall and precision at 0.75~\AA{}. Molecules are sorted by the recall difference and retain that ordering in both panels. Positive values favor Qwen, negative values favor LoQI. These descriptive differences are not per-molecule confidence intervals. '+common, 'HTML §7; drug_molecule_metrics.csv', True)

    fig, ax = plt.subplots(figsize=(8, 7.5))
    annotated_heatmap(ax, mat('pb_pass_rate')*100, names.loc[order].tolist(), labels)
    finish(fig, 'drug-validity', 'Chemical validity varies strongly across molecules and methods',
           r'PoseBusters pass percentages for each supplied drug molecule--method pool. The common color scale spans 0--100\%. Validity is reported separately from coverage and precision: coverage cannot be adjusted by multiplying it by these percentages because the passing conformers and recovered references need not coincide. '+common, 'HTML §7; drug_molecule_metrics.csv', True)

    fig, ax = plt.subplots(figsize=(8, 4.2))
    for y, key in enumerate(DK):
        row = drug[drug.label.eq(DRUG[key])]; a, b = row.mat_r.mean(), row.mat_p.mean()
        ax.plot([a, b], [y, y], c=STYLE[key][1]); ax.scatter(a, y, c=STYLE[key][1], marker='s', s=35); ax.scatter(b, y, edgecolors=STYLE[key][1], facecolors='white', s=35)
    method_axis(ax, DK); ax.set_xlabel('Mean matching distance (Å; lower is better)')
    ax.legend(handles=[Line2D([], [], marker='s', ls='', color='#444444', label='MAT-R: reference to generated'), Line2D([], [], marker='o', ls='', color='#444444', markerfacecolor='white', label='MAT-P: generated to reference')], loc='lower right', fontsize=7)
    finish(fig, 'drug-matching-distance', 'Continuous matching distances complement coverage',
           r'Molecule-mean MAT-R and MAT-P for the six selected generators. MAT-R averages nearest generated-conformer distances from observed references; MAT-P averages nearest observed-reference distances from generated conformers. Lower is better. Segments connect distinct summaries rather than uncertainty bounds. '+common, 'HTML §7; drug_molecule_metrics.csv', True)

    refs = drug.groupby('smiles').num_true_confs
    assert refs.nunique().eq(1).all()
    counts = refs.first().sort_values(); assert counts.sum() == 2450
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.barh(range(23), counts, color='#668daa'); ax.set_yticks(range(23), names.loc[counts.index]); ax.invert_yaxis()
    for y, n in enumerate(counts): ax.text(n+8, y, str(int(n)), va='center', fontsize=7)
    ax.set_xlim(0, counts.max()*1.12); ax.set_xlabel('Number of stored experimental reference conformers'); ax.grid(axis='x', alpha=.15); ax.set_axisbelow(True)
    finish(fig, 'drug-reference-counts', 'The experimental reference collection is unevenly sampled',
           r'Numbers of stored experimental reference conformers for the 23 drug molecules (2,450 in total). Reference counts are not counts of unique conformational states: redundant or very similar geometries may recur. Reported coverage and precision aggregates average molecules equally, so reference-rich molecules do not receive extra aggregate weight.', 'HTML §§1,7; drug_molecule_metrics.csv', True)

    indices = np.random.default_rng(SEED).integers(0, 23, (NBOOT, 23))
    other = [k for k in DK if k != 'loqi_raw']; rows=[]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.3), sharey=True)
    for ax, metric, title in zip(axes, ['cov_r_075', 'cov_p_075'], ['Recall difference (percentage points)', 'Precision difference (percentage points)']):
        matrix = matrices[metric]
        for y, key in enumerate(other):
            delta, low, high = bootstrap((matrix[DRUG[key]] - matrix[DRUG['loqi_raw']]).to_numpy()*100, indices)
            ax.errorbar(delta, y, xerr=[[delta-low], [high-delta]], fmt=STYLE[key][2], color=STYLE[key][1], capsize=3)
            rows.append(dict(key=key, metric=metric, delta=delta, ci_low=low, ci_high=high))
        ax.axvline(0, color='#555555', ls='--', lw=.8); ax.set_xlabel(title); ax.grid(axis='x', alpha=.15)
    axes[0].set_yticks(range(len(other)), [label(k) for k in other]); axes[0].invert_yaxis()
    pd.DataFrame(rows).to_csv(OUT / 'drug-paired-intervals.csv', index=False)
    finish(fig, 'drug-paired-uncertainty', 'How uncertain are differences from LoQI?',
           r'Paired differences from LoQI in molecule-mean COV-R and COV-P at 0.75~\AA{}. Error bars are exploratory 95\% percentile intervals from 10,000 paired molecule bootstrap samples (23 molecules; seed 20260928). These newly rendered intervals can differ slightly from earlier exports because the resampling seed differs. Intervals are not adjusted for multiple comparisons; an interval crossing zero does not establish equivalence. '+common, 'HTML §7; drug_molecule_metrics.csv', True)

    cases = ['Imatinib', 'Actinonin', 'Bortezomib', 'Dexamethasone']
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.5))
    for ax, name in zip(axes.flat, cases):
        subset = drug[drug.name.str.casefold().eq(name.casefold())].set_index('label')
        assert len(subset) == 6
        for key, shift in [('loqi_raw', -.18), (QWEN, .18)]:
            row = subset.loc[DRUG[key]]
            bars = ax.bar(np.arange(3)+shift, row[['cov_r_075','cov_p_075','pb_pass_rate']].to_numpy(float)*100, width=.34, color=STYLE[key][1], label=label(key))
            ax.bar_label(bars, fmt='%.1f', fontsize=6, padding=2)
        ax.set_title(f'{name} · {int(subset.num_true_confs.iloc[0])} references'); ax.set_xticks(range(3), ['Recall', 'Precision', 'PB pass']); ax.set_ylim(0, 110); ax.set_ylabel('Percent'); ax.grid(axis='y', alpha=.15); ax.set_axisbelow(True)
    axes[0, 0].legend(fontsize=7)
    finish(fig, 'drug-case-profiles', 'Illustrative failures differ across constraints',
           r'Four post hoc illustrative molecules comparing LoQI and Qwen 1.7B FSQ: recall, precision, and PoseBusters pass rate. Coverage uses the strict 0.75~\AA{} cutoff. The cases illustrate opposing recall preferences and disagreements between geometric coverage and chemical validity; they are not a representative sample or a structural explanation. '+common, 'HTML §7; drug_molecule_metrics.csv', True)


def write_catalog():
    (OUT / 'catalog.json').write_text(json.dumps(CATALOG, indent=2) + '\n')
    text = [r'% Generated by build_candidate_figures.py; review-only figure options.',
            r'\clearpage', r'\section*{Additional figure candidates for review}',
            r'These optional figures expand the analyses discussed in the HTML report using corrected archived data. They are included for figure selection, not as a proposed final figure count. The four main figures remain unchanged. Core analyses retain the single 0.75~\AA{} recovery cutoff; drug analyses remain exploratory. Historical energy and random-subsampling analyses are excluded.',
            r'\begingroup', r'\renewcommand{\thefigure}{C\ifnum\value{figure}<10 0\fi\arabic{figure}}', r'\setcounter{figure}{0}']
    for item in CATALOG:
        text += [r'\clearpage', r'\begin{figure}[p]', r'\centering',
                 r'\includegraphics[width=\linewidth,height=0.72\textheight,keepaspectratio]{figures/candidates/' + item['stem'] + '.pdf}',
                 r'\caption{\textbf{' + item['title'] + '.} ' + item['caption'] + '}',
                 r'\label{fig:candidate-' + item['stem'] + '}', r'\end{figure}']
    text += [r'\clearpage', r'\endgroup', '']
    (ROOT / 'figure-candidates.tex').write_text('\n'.join(text))
    md = ['# Additional figure candidates', '',
          '18 optional figures derived from the HTML report’s questions and corrected archived measurements. All are included at the end of the Overleaf draft through `figure-candidates.tex`. Comment out that input to hide the entire review appendix; remove individual figure blocks to select candidates. The four main figures are unchanged.', '',
          'Gallery: [additional-figures.pdf](figures/candidates/additional-figures.pdf). Every candidate also has PDF, SVG, and PNG artwork.', '',
          '| ID | Candidate | Report connection |', '| --- | --- | --- |']
    for item in CATALOG:
        md.append(f"| {item['id']} | [{item['title']}](figures/candidates/{item['stem']}.pdf) | {item['source']} |")
    md += ['', '## Reproduction and interpretation', '',
           '`python jcim-overleaf/build_candidate_figures.py` recreates the artwork and appendix using archived data only. Source hashes, cohort selection, and bootstrap settings are recorded in `figure-data/candidates/provenance.json`. Derived recovery matrices and confidence intervals are alongside the figures.', '',
           'C01–C10 describe core94 and the ten main generation pipelines, with ChEMBL3D-PB where applicable. Continuous RMSD and diversity summaries are conditional on measurable ensembles; binary recovery counts missing outputs as misses. C11–C18 describe the six selected generators on the 23-molecule drug set; coverage is not harmonized by PoseBusters filtering and failure handling. No historical energy curves, random-K results, additional recovery cutoffs, reference-cohort results, or training-recipe comparisons were introduced.', '',
           'These are alternatives as well as additions: C01/C02 emphasize uncertainty versus individual gains and losses; C03/C04 show individual outcomes versus pairwise complementarity; C11/C13/C18 show complete molecule profiles versus selected comparisons and examples. Keep the versions that serve the final story.', '']
    (ROOT / 'FIGURE_CANDIDATES.md').write_text('\n'.join(md))


def main():
    global GALLERY
    OUT.mkdir(exist_ok=True, parents=True)
    manifest = json.loads((DATA / 'candidates/provenance.json').read_text())
    assert hashlib.sha256((DATA / 'candidates/core-entries.csv').read_bytes()).hexdigest() == manifest['export_sha256']
    for src in manifest['sources']:
        assert hashlib.sha256((DATA / 'candidates' / src['path']).read_bytes()).hexdigest() == src['sha256']
    entries = pd.read_csv(DATA / 'candidates/core-entries.csv')
    tiers = pd.read_csv(DATA / 'recovery-table-tiers.csv')
    drug = pd.read_csv(DATA / 'drug_molecule_metrics.csv')
    with PdfPages(OUT / 'additional-figures.pdf') as gallery:
        GALLERY = gallery
        core_figures(entries, tiers)
        drug_figures(drug)
    assert len(CATALOG) == 18
    write_catalog()
    print(f'Rendered {len(CATALOG)} candidate figures and the review appendix.')


if __name__ == '__main__':
    main()
