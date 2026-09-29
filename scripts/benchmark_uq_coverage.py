#!/usr/bin/env python3
"""Repeated physical fixed-pilot coverage development benchmark.

Known finite Gaussian class and simulation bounds; not empirical calibration of
unknown experimental densities. All methods share each design and target. The
Bayesian baseline also receives an exactly matched prior-predictive check.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uncertainty import optimize_certificate,compress_nuisance_group
from fourier_splats.uq_baselines import gaussian_reference,prior_predictive_coverage,unregularized_reference
from fourier_splats.uq_physics import image_design,pose_jacobian,pose_remainder_bounds,pose_design_derivatives,perturbed_images,realify
from pilot_uncertainty import geometry,make_target

ROOT=Path(__file__).resolve().parents[1]


def run(seed,preferred,angle_deg,signals,replications):
    start=time.perf_counter();rng=np.random.default_rng(seed);n=96;p=40
    centers=rng.normal(size=(p//2,3))*2.2;centers[centers[:,2]<0]*=-1;sigma=1.
    pilot=rng.normal(size=p);pilot/=np.linalg.norm(pilot)
    B=.25;noise=.2;angle=np.deg2rad(angle_deg);shift=.05;box=24
    q,k,ctf=geometry(rng,n,preferred);ell=make_target(centers,sigma,preferred)
    a=image_design(k,ctf,centers,sigma,noise).reshape(-1,p)
    j=pose_jacobian(k,q,ctf,centers,sigma,pilot,box,angle,shift,noise)
    gamma,_,_=pose_remainder_bounds(k,q,ctf,centers,sigma,pilot,box,angle,shift,B,noise)
    gamma2,_,_=pose_remainder_bounds(k,q,ctf,centers,sigma,pilot,box,angle,shift,B,noise,structured=True)
    d=pose_design_derivatives(k,q,ctf,centers,sigma,box,angle,shift,noise)
    group=compress_nuisance_group(d.reshape(n,2*len(q[0]),-1))
    methods={};diagnostics={}
    for name,prior_multiplier,pose_scale in [('posterior',1,0),('posterior_wide',np.sqrt(p),0),('joint_gaussian',1,1/np.sqrt(5))]:
        prior_sd=B*prior_multiplier/np.sqrt(p)
        fit=gaussian_reference(a,ell,prior_sd,j,pose_scale)
        methods[name]=(fit['weights'],fit['posterior_half_width'])
        diagnostics[name]={'matched_linear_prior_predictive_coverage':prior_predictive_coverage(a,ell,prior_sd,fit,j,pose_scale),
                           'prior_sd':prior_sd,'pose_prior_sd':pose_scale}
        if name=='posterior':
            methods['diagonal_vi']=(fit['weights'],fit['diagonal_vi_half_width'])
            methods['ridge_sampling']=(fit['weights'],fit['measurement_noise_half_width'])
    ols=unregularized_reference(a,ell)
    if ols['identifiable']:methods['unregularized']=(ols['weights'],ols['half_width'])
    diagnostics['unregularized']={k:v for k,v in ols.items() if k!='weights'}
    for name,ji,eta,remainder,extra in [
        ('fixed_pose_bias_aware',np.zeros_like(j),0,0,[]),
        ('linearized_nuisance',j,1,0,[]),
        ('coarse_remainder',j,1,gamma,[]),
        ('structured_interaction',j,1,gamma2,[(group,B)]),
    ]:
        cert=optimize_certificate(a,ji,ell,B,eta,remainder,extra_nuisance_groups=extra,maxiter=200,rtol=1e-3)
        methods[name]=(cert.weights,cert.half_width)
        diagnostics[name]={'relative_gap':(cert.objective-cert.dual_lower_bound)/max(cert.objective,1e-15),
                           'converged':cert.converged,'iterations':cert.iterations,
                           'bias':cert.bias,'noise_sd':cert.noise_sd,'group_biases':cert.nuisance_group_biases}
    methods['no_data']=(np.zeros(a.shape[0]),B)
    names=list(methods);w=np.stack([methods[x][0] for x in names]);half=np.array([methods[x][1] for x in names])
    sd=np.linalg.norm(w,axis=1)
    # Factor the correlation matrix, not a covariance whose scales can differ
    # drastically when unregularized estimates amplify weak-view noise.
    normalized=w/np.maximum(sd[:,None],1e-300)
    correlation=normalized@normalized.T;ev,vec=np.linalg.eigh(correlation)
    root=sd[:,None]*vec*np.sqrt(np.maximum(ev,0))[None,:]
    rows=[]

    def evaluate(label,delta,u):
        expected=realify(perturbed_images(k,q,ctf,centers,sigma,pilot+delta,box,u,angle,shift,noise)).ravel()-a@pilot
        bias=w@expected-ell@delta
        analytic=np.ones(len(names))
        active=sd>1e-15
        analytic[active]=norm.cdf((half[active]-bias[active])/sd[active])-norm.cdf((-half[active]-bias[active])/sd[active])
        analytic[~active]=(np.abs(bias[~active])<=half[~active]+1e-12)
        errors=rng.normal(size=(replications,len(names)))@root.T+bias
        hits=(np.abs(errors)<=half+1e-12).sum(axis=0)
        for i,name in enumerate(names):
            rows.append({'truth_case':label,'method':name,'bias':float(bias[i]),'sd':float(sd[i]),
                         'half_width':float(half[i]),'width_fraction_of_no_data':float(half[i]/B),
                         'analytic_coverage':float(analytic[i]),'noise_replicates':replications,'covered':int(hits[i]),
                         'coefficient_bound_satisfied':bool(np.linalg.norm(delta)<=B+1e-12),
                         'pose_bound_satisfied':bool(np.max(np.linalg.norm(u,axis=1))<=1+1e-12)})

    for s in range(signals):
        delta=rng.normal(size=p);delta*=B*rng.uniform(.25,1)/np.linalg.norm(delta)
        u=rng.normal(size=(n,5));u*=rng.uniform(0,1,(n,1))/np.linalg.norm(u,axis=1,keepdims=True)
        evaluate(f'random_bounded_{s}',delta,u)
    for name,(wi,_) in methods.items():
        delta=a.T@wi-ell;delta*=B/max(np.linalg.norm(delta),1e-300)
        ui=np.einsum('nmq,nm->nq',j,wi.reshape(n,-1))
        ui/=np.maximum(np.linalg.norm(ui,axis=1,keepdims=True),1e-300)
        evaluate('adversarial_to_'+name,delta,ui)
    # Outside-theorem stress, never included in a claim of guaranteed coverage.
    wi=methods['posterior'][0];delta=a.T@wi-ell;delta*=2*B/np.linalg.norm(delta)
    ui=np.einsum('nmq,nm->nq',j,wi.reshape(n,-1));ui*=3/np.maximum(np.linalg.norm(ui,axis=1,keepdims=True),1e-300)
    evaluate('bounds_violated',delta,ui)
    return {'seed':seed,'geometry':'preferred' if preferred else 'uniform','angle_deg':angle_deg,
            'particles':n,'coefficients':p,'shift_pixels':shift,'density_radius':B,
            'seconds':time.perf_counter()-start,'diagnostics':diagnostics,'rows':rows}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seeds',type=int,default=12);p.add_argument('--signals',type=int,default=20)
    p.add_argument('--noise-replicates',type=int,default=2000);args=p.parse_args()
    target=ROOT/'results/uncertainty/development/coverage-benchmark.json'
    result={'stage':'development; fixed independent pilot, exact Gaussian class and known simulation bounds',
            'config':vars(args),'cases':[]}
    for seed in range(87100,87100+args.seeds):
        for preferred in [False,True]:
            for angle in [.5,2.]:
                case=run(seed,preferred,angle,args.signals,args.noise_replicates)
                result['cases'].append(case);target.write_text(json.dumps(result,indent=2)+'\n')
                summary={}
                for row in case['rows']:
                    if row['truth_case']!='bounds_violated':
                        summary[row['method']]=min(summary.get(row['method'],1),row['analytic_coverage'])
                print(seed,case['geometry'],angle,round(case['seconds'],2),summary,flush=True)
