"""Exact finite-pool recovery curves and the 80% recovery threshold.
Run after extract.py. Rejected candidates remain in the sampling population.
"""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
OUT=Path(__file__).resolve().parent
CFG=json.loads((OUT/'config.json').read_text())
def main():
    z=np.load(OUT/'candidate_rmsds.npz');dist=z['rmsd'];counts=z['n_candidates'];methods=z['methods'];ids=z['mol_ids']
    budgets=np.unique(np.rint(np.r_[np.geomspace(CFG['min_budget'],CFG['max_budget'],CFG['budget_grid_points']),CFG['explicit_budgets']]).astype(int))
    thresholds=np.array(CFG['rmsd_thresholds']);target=CFG['target_recovery'];n=len(ids)
    rng=np.random.default_rng(CFG['bootstrap_seed'])
    weights=rng.multinomial(n,np.full(n,1/n),size=CFG['bootstrap_replicates'])/n
    # Identical molecule resamples are used for every method, threshold and budget.
    lookup={}
    for N in np.unique(counts):
        if N==0:lookup[N]=np.zeros((len(budgets),1));continue
        K=np.minimum(budgets,N)
        lookup[N]=hypergeom.sf(0,int(N),np.arange(N+1)[None,:],K[:,None])
    curve_rows=[];frontier_rows=[];entry_rows=[];endpoint=[]
    probabilities=np.zeros((len(methods),len(thresholds),len(budgets),n))
    for mi,method in enumerate(methods):
        arrays=[dist[mi,i][np.isfinite(dist[mi,i])] for i in range(n)]
        for ti,t in enumerate(thresholds):
            P=np.stack([lookup[counts[mi,i]][:,np.searchsorted(arrays[i],t,side='right')] for i in range(n)],axis=1)
            probabilities[mi,ti]=P
            mean=P.mean(axis=1);boot=weights@P.T;low,high=np.quantile(boot,[.025,.975],axis=0)
            assert np.all(np.diff(mean)>=-1e-12)
            for j,k in enumerate(budgets):curve_rows.append(dict(method=str(method),budget=int(k),threshold=float(t),recovery=mean[j],ci_low=low[j],ci_high=high[j]))
            endpoint.append(dict(method=str(method),threshold=float(t),recovery=mean[-1]))
        # Exact inversion on the observed RMSD values; no arbitrary threshold grid.
        support=np.unique(np.concatenate(arrays));lo=np.zeros(len(budgets),int);hi=np.full(len(budgets),len(support)-1,int)
        def coverage(indices):
            ts=support[indices]
            return np.stack([lookup[counts[mi,i]][np.arange(len(budgets)),np.searchsorted(arrays[i],ts,side='right')] for i in range(n)]).mean(axis=0)
        attainable=coverage(hi)>=target
        while np.any(lo<hi):
            mid=(lo+hi)//2;above=coverage(mid)>=target
            hi=np.where(above,mid,hi);lo=np.where(above,lo,mid+1)
        frontier=support[lo];frontier[~attainable]=np.nan
        actual=coverage(lo);previous=coverage(np.maximum(lo-1,0))
        assert np.all(actual[attainable]>=target)
        assert np.all(previous[attainable & (lo>0)]<target)
        assert np.all(np.diff(frontier[np.isfinite(frontier)])<=1e-12)
        for j,k in enumerate(budgets):frontier_rows.append(dict(method=str(method),budget=int(k),rmsd_threshold=frontier[j],target_recovery=target,recovery_at_threshold=actual[j],attainable=bool(attainable[j])))
        print(method,'frontier',frontier[0],frontier[-1],flush=True)
    refs=z['reference_best_rmsd'];ref_threshold=np.sort(refs)[int(np.ceil(n*target))-1]
    pd.DataFrame(curve_rows).to_csv(OUT/'recovery_curves.csv',index=False)
    pd.DataFrame(frontier_rows).to_csv(OUT/'recovery80_frontier.csv',index=False)
    pd.DataFrame(endpoint).to_csv(OUT/'endpoints.csv',index=False)
    np.savez_compressed(OUT/'entry_probabilities.npz',probabilities=probabilities,methods=methods,mol_ids=ids,budgets=budgets,thresholds=thresholds)
    reference=dict(hit_rates={str(t):float(np.mean(refs<=t)) for t in thresholds},target_recovery=target,rmsd_threshold=float(ref_threshold),n_entries=n)
    (OUT/'reference.json').write_text(json.dumps(reference,indent=2)+'\n')
    print(reference,flush=True)
if __name__=='__main__':main()
