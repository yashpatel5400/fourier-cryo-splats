#!/usr/bin/env python3
"""Post-outcome classical comparators on immutable calibration arrays."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from fourier_splats.uq_view_risk import risk_baselines,paired_variance_decomposition
from fourier_splats.uq_moment_mc import rejection_count,projected_rejection_probability
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_view_risk.py','src/fourier_splats/uq_view_risk.py','tests/test_view_risk.py',
    'research/uncertainty/paired-power-v1/VIEW-RISK-BASELINES-PROTOCOL.md']
METHODS=['cvar_dkw','cvar_split','mean_only','unpaired_variance']


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():raise ValueError('Freeze source and protocol first')
    out=BASE/'bispectrum-view-risk-v1'
    if out.exists():raise ValueError('Preserve all outcomes')
    out.mkdir();result=dict(complete=False,cases=[],input_hashes={},sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    for ds in ['10028','10049']:
        p=BASE/f'bispectrum-view-variance-v1/{ds}/summary.json';r=json.loads(p.read_text());assert r['complete']
        ap=ROOT/r['array']['path'];assert sha(ap)==r['array']['sha256']
        for p in [p,ap]:result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        with np.load(ap) as f:
            for j,c in enumerate(r['cases']):
                for candidate in ['true','removed']:
                    center=c['calibration'][candidate]['center'];grouped=f[f'cal_{candidate}_grouped'][j]
                    baseline=risk_baselines(grouped,center,[1.,1.01,1.1,2.,5.])
                    decomposition=paired_variance_decomposition(grouped,center)
                    np.testing.assert_allclose(decomposition['variance_upper'],c['calibration'][candidate]['view_variance_upper'],atol=1e-15)
                    row=dict(dataset=ds,key=c['key'],candidate=candidate,baselines=baseline,decomposition=decomposition,projections=[])
                    for bounds in baseline['bounds']:
                        for method in METHODS:
                            for n in [1000,10000,100000]:
                                critical=rejection_count(n,bounds[method],.049)
                                old=next(p for p in c['projections'] if p['candidate']==candidate and p['method']=='view_variance' and p['kappa']==bounds['kappa'] and p['particles']==n)
                                row['projections'].append(dict(kappa=bounds['kappa'],method=method,particles=n,probability_bound=bounds[method],
                                    reject_at_count=critical,critical_minus_paired=critical-old['reject_at_count'],
                                    **{label:projected_rejection_probability(c['heldout'][label]['count'],r['heldout_views'],n,critical) for label in ['true','removed']}))
                    result['cases'].append(row)
                    print(ds,c['key'],candidate,'decomposition',json.dumps(decomposition),flush=True)
    result['complete']=True;(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print('COMPLETE',len(result['cases']),flush=True)


if __name__=='__main__':main()
