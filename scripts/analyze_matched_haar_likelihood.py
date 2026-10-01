#!/usr/bin/env python3
"""Matched numerical Haar mixture ratios; integration diagnostic, not a ceiling."""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
LEVELS=[4096,8192,16384,32768]
DELETIONS=[.1,.25,.5,1.]


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def realify(z):return np.ascontiguousarray(np.concatenate([z.real,z.imag],axis=1),dtype=float)


def mixture_diagnostics(scores,levels):
    """Stable prefix log means, importance ESS and largest normalized weight.

    Each quadrature point is exponentiated once. Accumulate disjoint chunks in
    log space so a late large kernel cannot underflow an earlier prefix.
    """
    x=np.asarray(scores,float)
    if x.ndim!=2 or not np.isfinite(x).all() or sorted(set(levels))!=list(levels) or levels[-1]>x.shape[1]:
        raise ValueError('Finite scores and increasing in-range prefix lengths required')
    logsum=np.full(len(x),-np.inf);logsum2=logsum.copy();maximum=logsum.copy()
    logs=[];esses=[];weights=[];begin=0
    for end in levels:
        block=x[:,begin:end];m=block.max(axis=1);e=np.exp(block-m[:,None])
        logsum=np.logaddexp(logsum,m+np.log(e.sum(axis=1)))
        logsum2=np.logaddexp(logsum2,2*m+np.log((e*e).sum(axis=1)))
        maximum=np.maximum(maximum,m)
        logs.append(logsum-np.log(end));esses.append(np.exp(2*logsum-logsum2));weights.append(np.exp(maximum-logsum))
        begin=end
    return np.array(logs),np.array(esses),np.array(weights)


def kernel_blocks(y0,region,m0,mr,deletion):
    """Dropped y-norm constants cancel in each within-image density ratio."""
    a=y0@m0.T;b=y0@mr.T;c=region@m0.T;d=region@mr.T
    q0=np.sum(m0*m0,axis=1);qr=np.sum(mr*mr,axis=1);cross=np.sum(m0*mr,axis=1)
    base=a-.5*q0
    delta=deletion
    change=-delta*b+delta*cross-.5*delta**2*qr
    return base,base+change,base-delta*c,base+change-delta*c+delta**2*d


def main():
    p=argparse.ArgumentParser();p.add_argument('--dataset',required=True,choices=['10028','10049','10076']);args=p.parse_args();ds=args.dataset
    sources=['scripts/analyze_matched_haar_likelihood.py','research/uncertainty/MATCHED-INFORMATION-LEDGER-PROTOCOL.md']
    for name in sources:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():raise ValueError('Commit source and protocol first')
    out=BASE/'matched-haar-information-v1'/ds
    if out.exists():raise ValueError('Preserve every attempted comparison')
    stage=BASE/'matched-information-ledger-v1'/ds;fisher=BASE/'candidate-fisher-score-v1'/ds
    a_summary=json.loads((stage/'summary.json').read_text());f_summary=json.loads((fisher/'summary.json').read_text())
    assert a_summary['complete'] and f_summary['complete']
    inputs={str(q.relative_to(ROOT)):sha(q) for q in [stage/'summary.json',stage/'replayed-images.npz',fisher/'summary.json',fisher/'arrays.npz']}
    assert inputs[str((stage/'replayed-images.npz').relative_to(ROOT))]==a_summary['replay']['sha256']
    assert inputs[str((fisher/'arrays.npz').relative_to(ROOT))]==f_summary['array']['sha256']
    out.mkdir(parents=True);started=time.perf_counter()
    result=dict(complete=False,dataset=ds,levels=LEVELS,deletions=DELETIONS,banks=2,views=8192,batch=64,
                numpy_version=np.__version__,python_version=platform.python_version(),
                thread_environment={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},
                source_hashes={p:sha(ROOT/p) for p in sources},inputs=inputs,
                git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),progress=[],cases=[],
                interpretation='Numerical quadrature under the fixed Haar/white-noise/one-CTF simulator; not a certified information ceiling or experimental likelihood.')
    def save():(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    save()
    with np.load(stage/'replayed-images.npz') as f:
        means=f['means'];noise=f['noise'];known_energy=f['region_energy']
    y0=realify(means[0]+noise);region=realify(means[1]);n=len(y0);assert n==8192
    with np.load(fisher/'arrays.npz') as f:
        train0=f['training_candidate_means'];trainr=f['training_region_means']
    ratios=np.empty((2,4,4,2,n));ess=np.empty((2,4,4,2,2,n));largest=np.empty_like(ess)
    for bank in range(2):
        m0=realify(train0[bank*32768:(bank+1)*32768]);mr=realify(trainr[bank*32768:(bank+1)*32768])
        q0=np.sum(m0*m0,axis=1);qr=np.sum(mr*mr,axis=1);cross=np.sum(m0*mr,axis=1)
        progress=dict(bank=bank,images_completed=0,complete=False,matmul_seconds=0.,quadrature_seconds=0.);result['progress'].append(progress);save()
        for begin in range(0,n,64):
            end=begin+64;yy=y0[begin:end];rr=region[begin:end]
            tick=time.perf_counter()
            a=yy@m0.T;b=yy@mr.T;c=rr@m0.T;d=rr@mr.T;base=a-.5*q0
            progress['matmul_seconds']+=time.perf_counter()-tick;tick=time.perf_counter()
            base_summary=mixture_diagnostics(base,LEVELS)
            for di,delta in enumerate(DELETIONS):
                change=-delta*b+delta*cross-.5*delta**2*qr
                alternatives=[mixture_diagnostics(base+change,LEVELS),mixture_diagnostics(base+change-delta*c+delta**2*d,LEVELS)]
                nulls=[base_summary,mixture_diagnostics(base-delta*c,LEVELS)]
                for state in [0,1]:
                    ratios[bank,:,di,state,begin:end]=alternatives[state][0]-nulls[state][0]
                    for model,info in enumerate([nulls[state],alternatives[state]]):
                        ess[bank,:,di,state,model,begin:end]=info[1]
                        largest[bank,:,di,state,model,begin:end]=info[2]
            progress['quadrature_seconds']+=time.perf_counter()-tick
            progress['images_completed']=end
            if end%512==0:
                save();print(ds,'bank',bank,'images',end,'seconds',round(time.perf_counter()-started,2),flush=True)
        progress['complete']=True;save()
    for bank in range(2):
        for li,level in enumerate(LEVELS):
            for di,delta in enumerate(DELETIONS):
                t0,t1=ratios[bank,li,di];gap=float(np.mean(t1-t0));v=float(np.var(t0,ddof=1))
                row=dict(bank=bank,quadrature_views=level,deletion=delta,null_logratio_mean=float(t0.mean()),
                         alternative_logratio_mean=float(t1.mean()),signed_mean_gap=gap,paired_gap_standard_error=float(np.std(t1-t0,ddof=1)/np.sqrt(n)),
                         null_variance=v,alternative_variance=float(np.var(t1,ddof=1)),squared_separation=gap*gap/v if v else None,
                         known_pose_twice_kl=float(delta**2*known_energy.mean()),diagnostics=[])
                for state in [0,1]:
                    for model in [0,1]:
                        ee=ess[bank,li,di,state,model];ww=largest[bank,li,di,state,model]
                        row['diagnostics'].append(dict(state=state,model=model,ess_min=float(ee.min()),
                            ess_quantiles={str(q):float(np.quantile(ee,q)) for q in [.01,.05,.5,.95]},
                            maximum_weight_quantiles={str(q):float(np.quantile(ww,q)) for q in [.5,.95,.99,1.]}))
                result['cases'].append(row)
    path=out/'likelihood-diagnostics.npz';np.savez_compressed(path,log_ratios=ratios,effective_samples=ess,maximum_weights=largest)
    result.update(complete=True,seconds=time.perf_counter()-started,arrays=dict(path=str(path.relative_to(ROOT)),sha256=sha(path),bytes=path.stat().st_size));save()
    print('COMPLETE',ds,result['seconds'],flush=True)


if __name__=='__main__':main()
