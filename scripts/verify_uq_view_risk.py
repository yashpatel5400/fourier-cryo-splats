#!/usr/bin/env python3
"""Independent scalar formula and minimal-critical-value replay for all risk baselines."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import binom
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def radius(x,width,delta):
    average=sum(x)/len(x);variance=sum((z-average)**2 for z in x)/(len(x)-1)
    return np.sqrt(2*variance*np.log(2/delta)/len(x))+7*width*np.log(2/delta)/(3*(len(x)-1))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'provenance/uncertainty/view-risk-independent-verification.json');args=parser.parse_args()
    out=args.output
    if out.exists():raise ValueError('Preserve previous verification')
    p=BASE/'bispectrum-view-risk-v1/summary.json';r=json.loads(p.read_text());assert r['complete']
    result=dict(complete=False,bound_checks=0,critical_value_checks=0,maximum_absolute_difference=0.,summary_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    for c in r['cases']:
        old=json.loads((BASE/f'bispectrum-view-variance-v1/{c["dataset"]}/summary.json').read_text());j=next(j for j,z in enumerate(old['cases']) if z['key']==c['key'])
        ap=ROOT/old['array']['path'];assert hashlib.sha256(ap.read_bytes()).hexdigest()==old['array']['sha256']
        with np.load(ap) as f:groups=f[f'cal_{c["candidate"]}_grouped'][j]
        x=groups.mean(axis=1);center=old['cases'][j]['calibration'][c['candidate']]['center'];m=len(x);mu=x.mean()
        rad=radius(x,1,.0005);lo=max(0,mu-rad);hi=min(1,mu+rad);closest=np.clip(.5,lo,hi);v=closest*(1-closest)
        joint=radius(x,1,.00025);jl=max(0,mu-joint);ju=min(1,mu+joint);sq=(x-center)**2
        vu=np.clip(sq.mean()+radius(sq,max(center**2,(1-center)**2),.0005)-max(0,jl-center,center-ju)**2,0,.25)
        eps=np.sqrt(np.log(2000)/(2*m));thresholds=np.unique(np.r_[0.,x,1.])
        first=np.sort(x[:m//2]);second=x[m//2:]
        for b in c['baselines']['bounds']:
            k=b['kappa'];dkw=min(1.,min(t+k*(np.maximum(0,x-t).mean()+eps*(1-t)) for t in thresholds))
            split_t=first[max(0,int(np.ceil((1-1/k)*len(first)))-1)];excess=np.maximum(0,second-split_t)
            expected=dict(cvar_dkw=dkw,cvar_split=min(1.,split_t+k*(excess.mean()+radius(excess,1-split_t,.001))),
                mean_only=min(1.,k*hi,hi+np.sqrt((k-1)*v)),unpaired_variance=min(1.,k*ju,ju+np.sqrt((k-1)*vu)))
            for name,value in expected.items():
                error=abs(b[name]-value);assert error<1e-12;result['maximum_absolute_difference']=max(result['maximum_absolute_difference'],error);result['bound_checks']+=1
        for p in c['projections']:
            k=p['reject_at_count'];n=p['particles'];prob=p['probability_bound']
            assert binom.sf(k-1,n,prob)<=.049+1e-13 and (k==0 or binom.sf(k-2,n,prob)>.049-1e-13);result['critical_value_checks']+=1
    result['complete']=True;out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
