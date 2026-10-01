#!/usr/bin/env python3
"""Predeclared classical comparators for the candidate-derived score study."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from fourier_splats.uq_view_risk import risk_baselines,paired_variance_decomposition
from fourier_splats.uq_moment_mc import rejection_count,projected_rejection_probability,binomial_interval
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/analyze_uq_candidate_scores.py','src/fourier_splats/uq_view_risk.py','research/uncertainty/paired-power-v1/CANDIDATE-SCORE-COMPARATORS.md']
METHODS=['cvar_dkw','cvar_split','mean_only','unpaired_variance']


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():raise ValueError('Freeze source and addendum first')
    out=BASE/'candidate-moment-score-v1/classical-comparators.json'
    if out.exists():raise ValueError('Preserve previous outcome')
    result=dict(complete=False,cases=[],input_hashes={},sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    for ds in ['10028','10049','10076']:
        sp=BASE/f'candidate-moment-score-v1/{ds}/summary.json';r=json.loads(sp.read_text());assert r['complete']
        ap=ROOT/r['array']['path'];assert sha(ap)==r['array']['sha256']
        for p in [sp,ap]:result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        with np.load(ap) as f:
            for j,c in enumerate(r['cases']):
                grouped=f['cal_grouped'][j];baseline=risk_baselines(grouped,.5,[1.,1.01,1.1,2.,5.])
                counts=np.asarray(c['counts']);events=f['heldout_scores'][j]>c['design']['threshold']
                row=dict(dataset=ds,key=c['key'],baselines=baseline,decomposition=paired_variance_decomposition(grouped,.5),projections=[],repeated_groups=[])
                for b in baseline['bounds']:
                    for method in METHODS:
                        for n in [1000,10000,100000]:
                            critical=rejection_count(n,b[method],.049);p=dict(kappa=b['kappa'],method=method,particles=n,probability_bound=b[method],reject_at_count=critical,outcomes=[])
                            for di,deletion in enumerate(r['deletions']):
                                for ai,a in enumerate(r['amplitudes']):p['outcomes'].append(dict(deletion=deletion,amplitude=a,**projected_rejection_probability(int(counts[di,ai]),r['heldout_views'],n,critical)))
                            row['projections'].append(p)
                        n=1024;critical=rejection_count(n,b[method],.049);p=dict(kappa=b['kappa'],method=method,particles=n,probability_bound=b[method],reject_at_count=critical,groups=128,outcomes=[])
                        for di,deletion in enumerate(r['deletions']):
                            for ai,a in enumerate(r['amplitudes']):
                                group_counts=events[di,ai].reshape(-1,n).sum(axis=1);rejected=int(np.sum(group_counts>=critical))
                                p['outcomes'].append(dict(deletion=deletion,amplitude=a,event_counts=group_counts.tolist(),rejected=rejected,fraction=rejected/128,pointwise_95_interval=binomial_interval(rejected,128)))
                        row['repeated_groups'].append(p)
                result['cases'].append(row);print(ds,c['key'],'classical comparators complete',flush=True)
    result['complete']=True;out.write_text(json.dumps(result,indent=2)+'\n');print('COMPLETE',len(result['cases']),flush=True)


if __name__=='__main__':main()
