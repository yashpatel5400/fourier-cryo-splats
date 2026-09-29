#!/usr/bin/env python3
"""Density-energy sensitivity certificates on real acquisition geometries.

An experimental pilot selects the density and its scale. Coverage labels come
from explicitly bounded Gaussian-dictionary simulations, not real-map truth.
No reference map chooses the radius. This is still a development experiment:
its class is an assumption and does not cover representation error.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_physics import (gaussian_pair_gram,density_energy_coordinates,
    gaussian_average_functionals,image_design,pose_jacobian,pose_design_derivatives,
    pose_remainder_bounds,perturbed_images,realify)
from fourier_splats.uncertainty import optimize_certificate,optimize_certificate_pdhg,compress_nuisance_group
from fourier_splats.uq_baselines import gaussian_reference

ROOT=Path(__file__).resolve().parents[1]


def run(args,dataset,n,angle_deg):
    started=time.perf_counter();rng=np.random.default_rng(args.seed)
    pilot_file=ROOT/'results/uncertainty/development/representation'/dataset/f'real_particles-spacing-{args.spacing}.npz'
    checkpoint=np.load(pilot_file);centers=checkpoint['centers'];sigma=float(checkpoint['sigma'])
    pilot=checkpoint['coefficients'].copy();gram=gaussian_pair_gram(centers,sigma)
    pilot/=np.sqrt(pilot@gram@pilot)  # exactly unit density energy by definition
    transform,energy_record=density_energy_coordinates(centers,sigma)
    operator_norm=1/np.sqrt(energy_record['smallest_retained_eigenvalue'])
    B=args.radius;angle=np.deg2rad(angle_deg);shift=args.shift if angle_deg else 0.
    g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=n,seed=args.seed)
    n=len(g['indices']);p=transform.shape[1];m=2*len(g['q'][0])
    raw=image_design(g['k'],g['ctf'],centers,sigma)
    pilot_images=raw@pilot
    noise=np.sqrt(np.mean(pilot_images**2)/args.snr)
    a=raw.reshape(-1,len(pilot))@transform/noise
    j=pose_jacobian(g['k'],g['q'],g['ctf'],centers,sigma,pilot,64,angle,shift,noise)
    # Transform of the density-pose interaction. Work in small particle blocks
    # so the uncompressed coefficient-by-pose tensor does not fill memory.
    blocks=[]
    if angle_deg:
        for start in range(0,n,8):
            sl=slice(start,start+8)
            d=pose_design_derivatives(g['k'][sl],g['q'][sl],g['ctf'][sl],centers,sigma,64,angle,shift,noise)
            energy_d=np.einsum('nmpq,pr->nmrq',d,transform,optimize=True)
            blocks.append(compress_nuisance_group(energy_d.reshape(len(d),m,-1)))
        group=np.concatenate(blocks);extra=[(group,B)]
    else:extra=[]
    # This transform norm bound is safe but can be loose for a dense dictionary.
    gamma,_,_=pose_remainder_bounds(g['k'],g['q'],g['ctf'],centers,sigma,pilot,64,angle,shift,B*operator_norm,noise,structured=True)
    points=np.array([[0,0,0],[.08,0,0],[0,0,.08],[0,0,-.08]])
    local=gaussian_average_functionals(points,centers,sigma,args.width)
    targets=[('center',local[0]),('x_offset',local[1]),('z_contrast',local[2]-local[3])]
    result={'stage':'development; real geometry and experimental pilot; assumed finite density-energy class',
            'dataset':dataset,'particles':n,'angle_degrees':angle_deg,'shift_pixels':shift,
            'snr_power':args.snr,'simulated_noise_std':float(noise),'density_energy_radius':B,
            'pilot_energy_norm':float(np.sqrt(pilot@gram@pilot)),'energy_coordinates':energy_record,
            'remainder_transform_norm':float(operator_norm),'frequency_radius_bins':5,
            'local_average_width_A':args.width*g['field_A'],'setup_seconds':time.perf_counter()-started,'targets':[]}
    for label,raw_ell in targets[:args.targets]:
        ell=transform.T@raw_ell;no_data=B*np.linalg.norm(ell);methods={};diagnostics={}
        for solver in args.solvers.split(','):
            call=optimize_certificate if solver=='irls' else optimize_certificate_pdhg
            t=time.perf_counter()
            fit=call(a,j,ell,B,1,gamma,extra_nuisance_groups=extra,
                     maxiter=args.maxiter if solver=='pdhg' else args.irls_maxiter,rtol=args.tolerance)
            methods[solver]=(fit.weights,fit.half_width)
            diagnostics[solver]={'half_width':fit.half_width,'width_fraction_of_no_data':fit.half_width/no_data,
                                 'noise_sd':fit.noise_sd,'density_bias':fit.reconstruction_bias,
                                 'nuisance_group_biases':fit.nuisance_group_biases,'remainder_bias':fit.remainder_bias,
                                 'relative_gap':(fit.objective-fit.dual_lower_bound)/max(fit.objective,1e-300),
                                 'converged':bool(fit.converged),'iterations':fit.iterations,'seconds':time.perf_counter()-t}
        posterior=gaussian_reference(a,ell,B/np.sqrt(p))
        methods['gaussian_posterior']=(posterior['weights'],posterior['posterior_half_width'])
        names=list(methods);w=np.stack([methods[x][0] for x in names]);widths=np.array([methods[x][1] for x in names]);sd=np.linalg.norm(w,axis=1)
        truths=[]
        for s in range(args.signals):
            delta=rng.normal(size=p);delta*=B/np.linalg.norm(delta)
            u=rng.normal(size=(n,5));u/=np.linalg.norm(u,axis=1,keepdims=True)
            truths.append((f'random_boundary_{s}',delta,u))
        for name,wi in zip(names,w):
            delta=a.T@wi-ell;delta*=B/max(np.linalg.norm(delta),1e-300)
            u=np.einsum('nmq,nm->nq',j,wi.reshape(n,m));u/=np.maximum(np.linalg.norm(u,axis=1,keepdims=True),1e-300)
            truths.append(('adversarial_to_'+name,delta,u))
        rows=[]
        for case,delta,u in truths:
            expected=realify(perturbed_images(g['k'],g['q'],g['ctf'],centers,sigma,pilot+transform@delta,64,u,angle,shift,noise)).ravel()-pilot_images.ravel()/noise
            bias=w@expected-ell@delta
            active=sd>1e-15;coverage=(np.abs(bias)<=widths+1e-12).astype(float)
            coverage[active]=norm.cdf((widths[active]-bias[active])/sd[active])-norm.cdf((-widths[active]-bias[active])/sd[active])
            for idx,name in enumerate(names):
                rows.append({'truth':case,'method':name,'bias':float(bias[idx]),'sd':float(sd[idx]),
                             'width_fraction_of_no_data':float(widths[idx]/no_data),'analytic_coverage':float(coverage[idx])})
        result['targets'].append({'target':label,'pilot_target':float(raw_ell@pilot),'no_data_half_width':float(no_data),
                                  'diagnostics':diagnostics,'coverage_rows':rows})
        print(dataset,n,angle_deg,label,diagnostics,flush=True)
    result['seconds']=time.perf_counter()-started
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--datasets',default='10028,10049,10076');parser.add_argument('--particles',default='128,512')
    parser.add_argument('--angles',default='0,0.5');parser.add_argument('--spacing',type=float,default=2.)
    parser.add_argument('--radius',type=float,default=.25);parser.add_argument('--snr',type=float,default=.1)
    parser.add_argument('--shift',type=float,default=.05);parser.add_argument('--width',type=float,default=.03)
    parser.add_argument('--solvers',default='pdhg');parser.add_argument('--maxiter',type=int,default=3000)
    parser.add_argument('--irls-maxiter',type=int,default=150);parser.add_argument('--tolerance',type=float,default=.002)
    parser.add_argument('--targets',type=int,default=3);parser.add_argument('--signals',type=int,default=10)
    parser.add_argument('--seed',type=int,default=609301);parser.add_argument('--output',default='real-geometry-energy.json')
    args=parser.parse_args();output=ROOT/'results/uncertainty/development'/args.output
    result={'config':vars(args),'cases':[]}
    for dataset in args.datasets.split(','):
        for n in map(int,args.particles.split(',')):
            for angle in map(float,args.angles.split(',')):
                result['cases'].append(run(args,dataset,n,angle))
                output.write_text(json.dumps(result,indent=2)+'\n')
