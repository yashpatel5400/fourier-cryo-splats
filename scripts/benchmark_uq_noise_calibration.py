#!/usr/bin/env python3
"""Noise-scale calibration controls for fixed ambient certificate weights.

Integrates inference Gaussian noise analytically. Calibration draws, not voxels,
are the independent replicates. Correlated calibration is a violation control.
"""
import json
from pathlib import Path
import numpy as np
from scipy.special import ndtr
from scipy.stats import norm
from fourier_splats.uq_noise import gaussian_scale_upper
from fourier_splats.uncertainty import bias_aware_half_width
ROOT=Path(__file__).resolve().parents[1]


def half_width(sd,bias,alpha):
    sd=np.asarray(sd);lo=np.full(sd.shape,bias);hi=bias+norm.ppf(1-alpha/2)*sd
    for _ in range(55):
        mid=(lo+hi)/2;coverage=ndtr((mid-bias)/sd)-ndtr((-mid-bias)/sd)
        less=coverage<1-alpha;lo=np.where(less,mid,lo);hi=np.where(less,hi,mid)
    return (lo+hi)/2


def run():
    rng=np.random.default_rng(609330);draws=20000;df=128;beta=.01;alpha=.04
    template=json.loads((ROOT/'results/uncertainty/development/ambient-pose.json').read_text())
    scenarios=[('independent_centered',0.,0.),('independent_signal_contaminated',.5,0.),('correlated_violation',0.,.5)]
    rows=[]
    for label,mean,correlation in scenarios:
        z=np.sqrt(1-correlation)*rng.normal(size=(draws,df))+np.sqrt(correlation)*rng.normal(size=(draws,1))+mean
        upper=gaussian_scale_upper(z,beta);plug=np.sqrt(np.mean(z*z,axis=1))
        for case in template['cases']:
            for angle in case['angles']:
                fit=angle['fits']['bounded_pose_ambient'];sd=fit['noise_sd']
                b=fit['density_bias']+sum(fit['nuisance_group_biases'])+fit['remainder_bias']
                methods={'known_sigma':np.full(draws,bias_aware_half_width(sd,b)),
                         'plug_in_sigma':half_width(plug*sd,b,.05),'independent_upper_sigma':half_width(upper*sd,b,alpha)}
                for name,q in methods.items():
                    if name!='known_sigma':
                        sample=3
                        np.testing.assert_allclose(q[sample],bias_aware_half_width((upper if name=='independent_upper_sigma' else plug)[sample]*sd,b,alpha if name=='independent_upper_sigma' else .05),rtol=1e-10)
                    for bias_fraction in [0.,1.]:
                        actual_bias=bias_fraction*b
                        cov=ndtr((q-actual_bias)/sd)-ndtr((-q-actual_bias)/sd)
                        rows.append({'scenario':label,'dataset':case['dataset'],'angle_deg':angle['angle_deg'],'method':name,
                                     'bias_fraction':bias_fraction,'coverage':float(cov.mean()),'calibration_mc_se':float(cov.std(ddof=1)/np.sqrt(draws)),
                                     'median_width_relative_known':float(np.median(q)/fit['half_width']),
                                     'scale_upper_failure_rate':float(np.mean(upper<1))})
        print(label,'scale failure',float(np.mean(upper<1)),flush=True)
    result={'stage':'development scalar noise calibration; fixed previously chosen ambient weights, not an estimated full covariance',
            'config':{'seed':609330,'calibration_draws':draws,'independent_coordinates':df,'alpha':alpha,'beta':beta,
                      'guaranteed_coverage_if_assumptions_hold':(1-alpha)*(1-beta)},'cases':rows}
    (ROOT/'results/uncertainty/development/noise-scale-calibration.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':run()
