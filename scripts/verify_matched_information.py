#!/usr/bin/env python3
"""Independent complex-moment and direct-distance checks of frozen diagnostics."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from fourier_splats.uq_bispectrum import MomentContrast
from verify_uq_bispectrum import IndependentComplexMoments

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results/uncertainty/development'


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    out = BASE / 'matched-information-verification-v1'
    if out.exists():
        raise ValueError('Preserve every verification attempt')
    out.mkdir()
    result = dict(complete=False, source_hash=sha(Path(__file__)), cases=[])
    for ds in ['10028', '10049', '10076']:
        ap = BASE / 'matched-information-ledger-v1' / ds
        bp = BASE / 'matched-haar-information-v1' / ds
        fp = BASE / 'candidate-fisher-score-v1' / ds
        aj = json.loads((ap/'summary.json').read_text())
        bj = json.loads((bp/'summary.json').read_text())
        fj = json.loads((fp/'summary.json').read_text())
        assert aj['complete'] and bj['complete'] and fj['complete']
        assert sha(ap/'replayed-images.npz') == aj['replay']['sha256']
        assert sha(bp/'likelihood-diagnostics.npz') == bj['arrays']['sha256']
        assert sha(fp/'arrays.npz') == fj['array']['sha256']
        with np.load(ap/'replayed-images.npz') as f:
            means, noise, energy = f['means'], f['noise'], f['region_energy']
        with np.load(fp/'arrays.npz') as f:
            train0, trainr, triads = f['training_candidate_means'], f['training_region_means'], f['triads']
            directions = {c['key']:f[c['key']+'_direction'] for c in fj['cases']}
        with np.load(bp/'likelihood-diagnostics.npz') as f:
            lr, ess, weight = f['log_ratios'], f['effective_samples'], f['maximum_weights']
        np.testing.assert_allclose(energy, (means[1].real**2+means[1].imag**2).sum(1), rtol=3e-15)
        row = dict(dataset=ds, input_hashes={str(p.relative_to(ROOT)):sha(p) for p in
            [ap/'summary.json',bp/'summary.json',fp/'summary.json',ap/'replayed-images.npz',bp/'likelihood-diagnostics.npz',fp/'arrays.npz']},
            variance_checks=[], direct_likelihood_checks=[], summary_checks=0)
        for case in aj['cases']:
            key=case['key']; basekey=key.removesuffix('_power_component').removesuffix('_bispectrum_component')
            w=directions[basekey].copy(); n=means.shape[-1]
            selected=triads if basekey.startswith('power_bispectrum') else np.empty((0,3),int)
            if key.endswith('_power_component'): w[n:]=0
            if key.endswith('_bispectrum_component'): w[:n]=0
            assert abs(w[:n].sum()-case['power_weight_sum']) < 1e-13
            model=MomentContrast(n,selected,w)
            enumerator=IndependentComplexMoments(n,selected,w)
            for index in [0,4095,8191]:
                for variance in [.9,1.,1.1]:
                    direct=enumerator.variance(means[0,index],variance)
                    analytic=model.mean_variance(means[0,index],noise_variance=variance)[1][0]
                    np.testing.assert_allclose(direct,analytic,rtol=2e-10,atol=2e-10)
                    row['variance_checks'].append(dict(key=key,index=index,noise_variance=variance,
                        direct=direct,analytic=float(analytic),error=float(abs(direct-analytic))))
            # Full null-variance decomposition recomputed from conditional moments.
            mm,vv=model.mean_variance(means[0])
            np.testing.assert_allclose([mm.mean(),vv.mean(),mm.var(),vv.mean()+mm.var()],
                [case[k] for k in ['null_mean','noise_variance_mean','view_mean_variance','null_total_variance']],rtol=1e-13,atol=1e-13)
            for dc in case['deletions']:
                np.testing.assert_allclose(dc['known_pose_twice_kl'],dc['deletion']**2*energy.mean(),atol=1e-15)
                np.testing.assert_allclose(dc['squared_separation'],dc['signed_mean_gap']**2/case['null_total_variance'],atol=1e-15)
            row['summary_checks']+=1
        # All stored summary statistics are replayed from the raw likelihood arrays.
        for case in bj['cases']:
            bi=case['bank']; li=bj['levels'].index(case['quadrature_views']); di=bj['deletions'].index(case['deletion'])
            t0,t1=lr[bi,li,di]; gap=np.mean(t1-t0); vv=np.var(t0,ddof=1)
            np.testing.assert_allclose([t0.mean(),t1.mean(),gap,vv,gap*gap/vv],
                [case[k] for k in ['null_logratio_mean','alternative_logratio_mean','signed_mean_gap','null_variance','squared_separation']],atol=1e-14,rtol=1e-14)
            e=ess[bi,li,di]; ww=weight[bi,li,di]
            assert np.isfinite(e).all() and np.isfinite(ww).all()
            assert e.min() >= 1-1e-10 and e.max() <= case['quadrature_views']+1e-7
            assert ww.min() >= 1/case['quadrature_views']-1e-12 and ww.max() <= 1+1e-10
            row['summary_checks']+=1
        # Direct complex Euclidean residuals, no dropped norms or four-dot-product identity.
        for bank in [0,1]:
            start=bank*32768; m0=train0[start:start+32768]; mr=trainr[start:start+32768]
            for index in [0,4095,8191]:
                y0=means[0,index]+noise[index]
                for di,delta in enumerate(bj['deletions']):
                    for state in [0,1]:
                        y=y0-state*delta*means[1,index]
                        kernels=[-.5*np.sum(abs(y-mm)**2,axis=1) for mm in [m0,m0-delta*mr]]
                        for li,level in enumerate(bj['levels']):
                            logs=[logsumexp(k[:level]) for k in kernels]
                            actual=logs[1]-logs[0]; error=abs(actual-lr[bank,li,di,state,index])
                            if error>1e-9: raise ArithmeticError('Direct Gaussian mixture replay failed')
                            diagnostics=[]
                            for model in [0,1]:
                                probs=np.exp(kernels[model][:level]-logs[model]); e=1/(probs@probs); w=probs.max()
                                np.testing.assert_allclose(e,ess[bank,li,di,state,model,index],rtol=2e-10,atol=2e-10)
                                np.testing.assert_allclose(w,weight[bank,li,di,state,model,index],rtol=2e-10,atol=2e-10)
                                diagnostics.append(dict(ess=float(e),maximum_weight=float(w)))
                            row['direct_likelihood_checks'].append(dict(bank=bank,index=index,deletion=delta,state=state,level=level,error=float(error),diagnostics=diagnostics))
        row['maximum_variance_error']=max(c['error'] for c in row['variance_checks'])
        row['maximum_direct_likelihood_error']=max(c['error'] for c in row['direct_likelihood_checks'])
        result['cases'].append(row)
        (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
        print(ds,len(row['variance_checks']),len(row['direct_likelihood_checks']),row['maximum_variance_error'],row['maximum_direct_likelihood_error'],flush=True)
    result['complete']=True
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
