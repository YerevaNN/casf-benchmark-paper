"""Verify only the presented landscape inputs, annotations, and reference counts."""
from pathlib import Path
import hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent;PAPER=HERE.parents[1]
manifest=json.loads((HERE/'provenance.json').read_text())
for section in ['inputs','outputs']:
    for p,h in manifest[section].items():assert hashlib.sha256((PAPER/p).read_bytes()).hexdigest()==h,p
meta=json.loads((HERE/'plot-data.json').read_text());land=json.loads((HERE/'landscape-data.json').read_text())
entries=json.loads((PAPER/'plinder22-data/entries.json').read_text());lookup={(r['method'],r['smiles']):r['key']for r in entries}
original=np.load(PAPER/'plinder22-data/nearest-distances.npz',allow_pickle=False)
def labels(d,order):
    centers=[];out=np.empty(len(d),dtype=int)
    for i in order:
        if centers:
            j=int(np.argmin(d[i,centers]))
            if d[i,centers[j]]<.5:out[i]=j;continue
        out[i]=len(centers);centers.append(i)
    return out
with np.load(HERE/'reference_distances.npz',allow_pickle=False)as z:
    mats={str(s):z[f'm{i}']for i,s in enumerate(z['smiles'])}
assert set(meta['names'].values())=={'Staurosporine','Actinonin'}
for s in meta['order']:
    d=mats[s];assert np.isfinite(d).all()and np.allclose(d,d.T)and np.allclose(np.diag(d),0)
    rng=np.random.default_rng(42+meta['molecule_ordinals'][s]);perms=[np.arange(len(d))]+[rng.permutation(len(d))for _ in range(20)]
    for pi,perm in enumerate(perms):
        counts=np.bincount(labels(d,perm));r=next(r for r in meta['clusters']if r['name']==meta['names'][s]and r['permutation']==pi)
        assert r['n_references']==len(d)and r['n_clusters']==len(counts)and r['largest_cluster']==counts.max()
for m in land['molecules']:
    s=m['smiles'];d=m['dedup'];assert len(m['pts']['s'])==len(d['pts']['x'])==len(d['pts']['y'])
    assert len(d['ref_idx'])==len(mats[s])
    for pagekey,method in [('loqi','loqi_druglike'),('flowr','flowr_druglike'),('qwen_v2_c17','qwen_1p7b_v2_coordtoken_step41054')]:
        nn=np.array(d['nn_ref'][pagekey]);ref=original[lookup[method,s]+'_reference'];gen=original[lookup[method,s]+'_generated']
        assert np.isfinite(nn).all()and len(nn)==len(gen)and len(ref)==len(d['ref_idx'])
        np.testing.assert_allclose((nn<.75).mean(),(gen<=.75).mean(),atol=1e-12,rtol=0)
print('Verified two landscape examples, six source panels, canonical metric annotations, 42 reference clusterings, and source/artwork hashes.')
