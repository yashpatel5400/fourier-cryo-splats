#!/usr/bin/env python3
"""Independent post hoc checks; preserve raw author-code outputs and diagnostics."""
import hashlib
import json
from pathlib import Path
import numpy as np
from replay_population_baseline import ROOT, AUTHOR, OUT, load


def main():
    destination = OUT / 'independent-check.json'
    if destination.exists():
        raise ValueError('Preserve previous check')
    summary = json.loads((OUT/'summary.json').read_text())
    assert summary['complete'] and len(summary['cases']) == 23
    rows = []
    for row in summary['cases']:
        p = AUTHOR / row['input']
        assert hashlib.sha256(p.read_bytes()).hexdigest() == summary['inputs'][row['input']]['sha256']
        ell = load(p)['log_likelihoods'] if row['kind']=='synthetic' else np.load(p,allow_pickle=False)
        ell = np.asarray(ell, dtype=np.longdouble)
        # tanh log-odds and bisection use different arithmetic from exp/brentq.
        z = np.tanh((ell[:,0]-ell[:,1])/2)
        lo, hi = np.longdouble(0), np.longdouble(1)
        for _ in range(80):
            mid = (lo+hi)/2
            derivative = np.mean(2*z/(1+(2*mid-1)*z))
            if derivative>0: lo=mid
            else: hi=mid
        t = float((lo+hi)/2)
        w = np.asarray(row['author_weights'],float)
        normalization_error = float(w.sum()-1)
        w = w/w.sum()
        ll = np.asarray(ell,float)
        q = np.exp(ll-ll.max(axis=1,keepdims=True))
        mean_loss = float(np.mean(np.log(q@np.array([t,1-t]))-np.log(q@w)))
        gap = float(np.max(np.mean(q/(q@w)[:,None],axis=0))-1)
        difference = abs(t-row['independent_mle'][0])
        if difference>1e-10 or mean_loss < -1e-12:
            raise ArithmeticError('Independent optimizer check failed')
        rows.append(dict(input=row['input'],root_difference=difference,
                         author_weight_sum_error=normalization_error,
                         normalized_author_weights=w.tolist(),
                         normalized_weight_difference=float(np.max(abs(w-[t,1-t]))),
                         normalized_mean_objective_loss=mean_loss,
                         normalized_simplex_gap=gap))
    equal = []
    for i,r in enumerate(summary['cases']):
        for s in summary['cases'][i+1:]:
            if summary['inputs'][r['input']]['sha256']==summary['inputs'][s['input']]['sha256']:
                equal.append([r['input'],s['input']])
    result=dict(complete=True,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                summary_sha256=hashlib.sha256((OUT/'summary.json').read_bytes()).hexdigest(),
                longdouble_precision_bits=int(np.finfo(np.longdouble).nmant),
                cases=rows,identical_likelihood_files=equal,
                note='Raw diagnostics use unnormalized returned JAX weights. Their tiny simplex drift can cause negative objective differences; use the normalized checks here for optimization comparisons. Apple longdouble precision is reported, not assumed extended.')
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['complete','longdouble_precision_bits','identical_likelihood_files']}))
    print('maximum root difference',max(x['root_difference'] for x in rows))


if __name__=='__main__':main()
