"""Validate selected SI against current main-text data and LaTeX inclusion order."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
PAPER=HERE.parent
manifest=json.loads((HERE/'provenance.json').read_text())
for section in ['inputs','outputs']:
    for name,digest in manifest[section].items():
        assert hashlib.sha256((PAPER/name).read_bytes()).hexdigest()==digest,name

changes=pd.read_csv(HERE/'recovery-gains-losses.csv')
assert len(changes)==20
np.testing.assert_array_equal(changes[['gained','lost','shared_hits','shared_misses']].sum(axis=1),changes.n)
assert changes.loc[(changes.cohort=='ref')&(changes.name=='Qwen 1.7B FSQ'),['gained','lost']].iloc[0].tolist()==[155,15]

counts=pd.read_csv(PAPER/'results-data/three-dataset-tables/table3.csv')
styles=json.loads((PAPER/'sampling-data/palette.json').read_text())
for cohort,path in [('core','figure-data/candidates/core-entries.csv'),('ref','results-data/ref-panel/per-entry.csv')]:
    source=pd.read_csv(PAPER/path)
    for key,(name,*_) in styles.items():
        for target in (['chembl_count'] if key=='chembl3d_gt_pb' else ['chembl_count','fixed']):
            method=key if key=='chembl3d_gt_pb' else key+'_'+target
            rows=source[source.method==method]
            assert len(rows)==(94 if cohort=='core' else 1044)
            expected=counts[(counts.cohort==cohort)&(counts.method==name)&(counts.budget==target)].iloc[0]
            np.testing.assert_allclose(rows.post_pb_confs.fillna(0).mean(),expected.mean_retained,atol=1e-10)
            tex=(PAPER/'tables/supporting-ensemble-counts.tex').read_text()
            line=next(line for line in tex.splitlines() if line.startswith(name+' &'))
            index=(0 if cohort=='core' else 2)+(0 if target=='chembl_count' else 1)
            assert line.split(' & ')[index+1].replace(r' \\','')==f'{expected.mean_retained:.1f}'

pl=pd.read_csv(HERE/'plinder-by-molecule.csv')
assert len(pl)==242 and pl.label.nunique()==11 and pl.smiles.nunique()==22
assert not pl.label.isin(['torsion_raw_druglike','torsion_minimized_22_druglike']).any()
q=pd.read_csv(HERE/'qwen-validity-source.csv')
assert len(q)==88 and q.model.nunique()==4
matched=pl[pl.label=='qwen_1p7b_v2_coordtoken_step41054'].merge(q[q.model=='fsq17'],on='smiles',validate='one_to_one',suffixes=('_main','_comparison'))
for column in ['n_generated','n_pb_pass']:
    np.testing.assert_array_equal(matched[column+'_main'],matched[column+'_comparison'])

# Expand includes in place to check the automatic numbering in the separate SI.
def expand(path):
    text=re.sub(r'(?<!\\)%[^\n]*','',path.read_text())
    def include(match):
        name=match[1].replace(r'\manuscriptroot','').strip()
        return expand(PAPER/(name if name.endswith('.tex') else name+'.tex'))
    return re.sub(r'\\input\{([^}]+)\}',include,text)
si=expand(PAPER/'supporting-information.tex')
expected_tables=['supporting-ensemble-counts','plinder23-screen','plinder23-statistics','plinder23-molecules','generators','energy-summary','selection-counts','qwen-training']
tables=re.findall(r'\\label\{tab:([^}]+)\}',si)
# Exact IDs are retained from the existing source files.
assert tables==expected_tables,tables
figures=re.findall(r'\\label\{fig:([^}]+)\}',si)
assert len(figures)==6,figures
assert figures[1:]==['qwen-molecule-effects','supporting-gains-losses','supporting-composition','supporting-plinder-coverage','supporting-plinder-validity'],figures
assert 'figures/candidates/' not in si
for forbidden in ['Mean median','tab:supporting-druglike','tab:size','tab:supporting-rmsd']:
    assert forbidden not in si
main=expand(PAPER/'acs-template.tex')
assert 'figures/candidates/' not in main and 'Additional figure candidates for review' not in main
assert expand(PAPER/'supporting-content.tex') in main
assert expand(PAPER/'supporting-content.tex') in si
for label in ['tab:'+x for x in tables]+['fig:'+x for x in figures]:
    assert main.count(r'\label{'+label+'}')==1,label
for source in [main,si]:
    for kind,number in re.findall(r'(?:Supporting )?(Table|Figure)s?\s*~?\s*S(\d+)',source):
        assert 1<=int(number)<=(8 if kind=='Table' else 6),(kind,number)
    for name in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',source):
        assert (PAPER/name).is_file(),name
print('Verified current source/artwork hashes, 20 paired recovery comparisons, retained counts, 242 PLINDER and 88 v2 validity records.')
print('Current SI: eight tables and six figures; no historical review gallery or removed tables included.')
print('SI table order:',tables)
print('SI figure order:',figures)

# Selection counts come from executed logs, not reconstructed exclusions.
audit=PAPER/'method-data/selection-audit'
for metadata in json.loads((audit/'sources.json').read_text()).values():
    path=audit/metadata['archive']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==metadata['sha256']
expected=[('casf_ligands_processed',285,4838),('casf_regen_failures',0,11),
          ('excluded_no_smiles_match',189,3566),('excluded_missing_topology_sdf',1,12),
          ('excluded_no_rotatable_bonds',1,13),('exact_matched_ligands',94,1236)]
tex=(PAPER/'selection-appendix.tex').read_text()
for key,core,ref in expected:
    for cohort,value in [('core',core),('ref',ref)]:
        log=(audit/f'{cohort}-matching-log.txt').read_text()
        assert re.search(r'^'+key+'='+str(value)+r'$',log,re.M),(cohort,key)
    assert f' & {core:,} & {ref:,} ' in tex,key
assert 285-0-189-1-1==94
assert 4838-11-3566-12-13==1236
selection=json.loads((PAPER/'results-data/ref-panel/selection.json').read_text())
assert selection['excluded_entries']==192 and selection['selected_entries']==1044
assert ' & 0 & 192 ' in tex and ' & 94 & 1,044 ' in tex
print('Verified selection table against archived execution logs and frozen ref selection.')

# Check training statements against completed-run records and upstream filter receipts.
audit=PAPER/'method-data/training-audit'
for metadata in json.loads((audit/'sources.json').read_text())['archived'].values():
    assert hashlib.sha256((audit/metadata['archive']).read_bytes()).hexdigest()==metadata['sha256']
filtered=json.loads((audit/'filter-summary.json').read_text())
removed=sum(filtered[k] for k in ['casf_similar','plinder_only_dropped','casf_ref_exact'])
assert filtered['rows']-removed==filtered['kept']==187398101
assert removed==2393927
per_source=pd.read_csv(audit/'filter-per-source.csv')
for key in ['rows','kept','casf_similar','plinder_only_dropped','casf_ref_exact']:
    assert int(per_source[key].sum())==filtered[key],key
for epoch in [1,2]:
    manifest=json.loads((audit/f'coordtoken-epoch{epoch}-manifest.json').read_text())
    assert manifest['data_manifest_sha256']==hashlib.sha256((audit/'coordtoken-input-manifest.json').read_bytes()).hexdigest()
    result=json.loads((audit/f'coordtoken-epoch{epoch}-result.json').read_text())
    assert result['status']=='complete' and result['steps']==366012
    assert result['examples_seen']==filtered['kept']
encoded=json.loads((audit/'fsq-encoding-manifest.json').read_text())
inputs=json.loads((audit/'coordtoken-input-manifest.json').read_text())
assert encoded['source']==inputs['sources']['train']['path']
release=json.loads((audit/'coordtoken-release.json').read_text())
assert encoded['checkpoint_sha256']==release['checkpoint_sha256']
tex=(PAPER/'training-appendix.tex').read_text()
configs=json.loads((PAPER/'results-data/qwen-size-representation/training-configurations.json').read_text())
generations=json.loads((PAPER/'results-data/qwen-size-representation/generation-metadata.json').read_text())
for key,label in [('binned06','0.6B binned'),('binned17','1.7B binned'),('fsq06','0.6B FSQ'),('fsq17','1.7B FSQ')]:
    final=json.loads((audit/f'{key}-final-update.json').read_text())
    exposure=final['exposure']
    assert exposure['completed_epochs']==1
    assert final['update']==configs[key]['config']['training']['optimizer_updates']
    row=' & '.join([label]+[f'{exposure[k]:,}' for k in ['documents','completed_updates','real_input_tokens','supervised_targets']])
    assert row in tex,row
    representation='binned' if key.startswith('binned') else 'fsq'
    dataset=json.loads((audit/f'{representation}-dataset.json').read_text())
    assert exposure['documents']==dataset['counts']['train']['accepted']
    for cohort,metadata in generations[key].items():
        shards=metadata.get('shards',[metadata])
        for shard in shards:
            assert shard['temperature']==shard['top_p']==1
            assert shard['samples_per_molecule']==1000
            assert shard['seed']==42+shard.get('shard_index',0)
            if representation=='fsq':
                assert shard['coordtoken_checkpoint_sha256']==encoded['checkpoint_sha256']
print('Verified training filter totals, source counts, tokenizer lineage, completed Qwen exposure and recorded generation settings.')
