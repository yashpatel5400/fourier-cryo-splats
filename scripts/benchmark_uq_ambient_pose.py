#!/usr/bin/env python3
"""Combined supported-voxel and bounded nonlinear pose development experiment."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import norm
from fourier_splats.uq_data import (particle_geometry,VoxelObservationOperator,SupportedVoxelOperator,
    voxel_pose_terms,local_weights,VoxelReference)
from fourier_splats.uq_physics import density_functionals
from fourier_splats.uq_subspace import matrix_free_certificate
from fourier_splats.uncertainty import optimize_certificate
from fourier_splats.uq_baselines import gaussian_reference_operator,prior_predictive_coverage

ROOT=Path(__file__).resolve().parents[1]
MAPS={'10028':'2660','10049':'6487','10076':'8434'}


def run(args,dataset):
    rng=np.random.default_rng(args.seed);start=time.perf_counter()
    g=particle_geometry(ROOT,dataset,'inference_half0',radius=args.frequency_radius,count=args.particles,seed=args.seed)
    box=args.box;grid=np.arange(-box//2,box//2)/box;z,y,x=np.meshgrid(grid,grid,grid,indexing='ij')
    xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1);mask=np.sum(xyz*xyz,axis=1)<=args.support**2
    ck=np.load(ROOT/'results/uncertainty/development/representation'/dataset/'real_particles-spacing-2.0.npz')
    pilot=density_functionals(xyz[mask],ck['centers'],float(ck['sigma']))@ck['coefficients'];pilot/=np.linalg.norm(pilot)
    raw=SupportedVoxelOperator(VoxelObservationOperator(g['k'],g['ctf'],box),mask)
    noise=np.sqrt(np.mean(raw.forward(pilot)**2)/args.snr)
    op=SupportedVoxelOperator(VoxelObservationOperator(g['k'],g['ctf'],box,noise),mask)
    a=op.as_linear_operator();n=op.op.n;m=2*op.op.q;B=args.radius
    reference=VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map',box=box).volume.ravel()[mask]
    reference/=np.linalg.norm(reference)
    if args.target=='center':
        ell=local_weights(box,[0,0,0],args.width*box).ravel()[mask]
    else:
        ell=(local_weights(box,[0,0,args.offset*box],args.width*box)-local_weights(box,[0,0,-args.offset*box],args.width*box)).ravel()[mask]
    no_data=B*np.linalg.norm(ell);baseline=matrix_free_certificate(op,ell,B,rtol=args.tolerance)
    result={'stage':'development; joint finite supported-voxel/pose assumptions, known simulated noise',
            'dataset':dataset,'config':vars(args),'particles':n,'active_voxels':len(pilot),'noise_std':float(noise),
            'reference_pilot_distance':float(np.linalg.norm(reference-pilot)),
            'no_data_half_width':float(no_data),'setup_seconds':time.perf_counter()-start,'angles':[]}
    for angle_deg in map(float,args.angles.split(',')):
        angle=np.deg2rad(angle_deg);shift=args.shift if angle_deg else 0.
        begin=time.perf_counter();terms=voxel_pose_terms(op,pilot,g['q'],angle,shift,B,order=2 if args.quadratic_pose else 1)
        setup_seconds=time.perf_counter()-begin;j=terms['jacobian'];group=terms['interaction_factor']
        gamma=.5*(terms['pilot_second_bound']+B*terms['operator_second_bound'])
        methods={'fixed_pose_ambient':(baseline.weights,baseline.half_width)};fits={}
        configurations=[('linearized_ambient',0,[]),('bounded_pose_ambient',gamma,[(group,B)])]
        if args.quadratic_pose:
            configurations.append(('quadratic_pose_ambient',terms['remainder'],[(group,B),
                (terms['quadratic_jacobian'].reshape(n,m,25),.5),(terms['quadratic_interaction_factor'],.5*B)]))
        for name,remainder,extra in configurations:
            begin=time.perf_counter();fit=optimize_certificate(a,j,ell,B,1,remainder,extra_nuisance_groups=extra,
                gram_diagonal=op.gram_diagonal,maxiter=args.maxiter,rtol=args.tolerance,cg_maxiter=500,cg_rtol=1e-7)
            methods[name]=(fit.weights,fit.half_width)
            fits[name]={'half_width':float(fit.half_width),'width_fraction_of_no_data':float(fit.half_width/no_data),
                        'noise_sd':fit.noise_sd,'density_bias':fit.reconstruction_bias,'nuisance_group_biases':fit.nuisance_group_biases,
                        'remainder_bias':fit.remainder_bias,'objective':fit.objective,'dual_lower_bound':fit.dual_lower_bound,
                        'relative_gap':float((fit.objective-fit.dual_lower_bound)/fit.objective),'converged':fit.converged,
                        'iterations':fit.iterations,'seconds':time.perf_counter()-begin,'history':fit.history}
            print(dataset,angle_deg,name,{k:v for k,v in fits[name].items() if k!='history'},flush=True)
        gaussian_checks=[]
        if args.gaussian_baselines:
            for factor in [.25,1.,4.]:
                prior_sd=factor*B/np.sqrt(len(pilot))
                for pose_sd in [0.,1/np.sqrt(5)]:
                    name=f'Gaussian_prior_energy_{factor}_pose_{pose_sd:.3f}'
                    reference_fit=gaussian_reference_operator(a,ell,prior_sd,j,pose_sd,gram_diagonal=op.gram_diagonal)
                    methods[name]=(reference_fit['weights'],reference_fit['posterior_half_width'])
                    gaussian_checks.append({'method':name,'density_prior_sd':prior_sd,'pose_prior_sd':pose_sd,
                        'matched_additive_linear_prior_predictive_coverage':prior_predictive_coverage(a,ell,prior_sd,reference_fit,j,pose_sd),
                        'relative_normal_residual':reference_fit['relative_normal_residual']})
        names=list(methods);weights=np.stack([methods[name][0] for name in names]);half=np.array([methods[name][1] for name in names]);sd=np.linalg.norm(weights,axis=1)
        cases=[]
        for idx in range(args.signals):
            delta=rng.normal(size=len(pilot));delta*=B/np.linalg.norm(delta)
            u=rng.normal(size=(n,5));u/=np.linalg.norm(u,axis=1,keepdims=True)
            cases.append((f'random_boundary_{idx}',delta,u))
        for name,w in zip(names,weights):
            delta=op.adjoint(w)-ell;delta*=B/max(np.linalg.norm(delta),1e-300)
            u=np.einsum('nmq,nm->nq',j,w.reshape(n,m));u/=np.maximum(np.linalg.norm(u,axis=1,keepdims=True),1e-300)
            cases.append(('adversarial_to_'+name,delta,u))
        coherent=np.zeros((n,5));coherent[:,0]=1
        cases.append(('independent_EMDB_reference_coherent_pose',reference-pilot,coherent))
        # Target-directed signals reveal power, rather than only coverage at random targets.
        for fraction in [.1,.25,.5,1.]:
            delta=fraction*B*ell/np.linalg.norm(ell)
            cases.append((f'target_direction_fraction_{fraction}',delta,np.zeros((n,5))))
        records=[]
        for label,delta,u in cases:
            r=Rotation.from_rotvec(angle*u[:,:3]).as_matrix();kp=np.einsum('nqi,nij->nqj',g['k'],r)
            perturbed=SupportedVoxelOperator(VoxelObservationOperator(kp,g['ctf'],box,noise),mask)
            real=perturbed.forward(pilot+delta).reshape(n,m);values=real[:,:m//2]+1j*real[:,m//2:]
            values*=np.exp(-2j*np.pi*shift*np.einsum('nqi,ni->nq',g['q'],u[:,3:])/box)
            expected=np.concatenate([values.real,values.imag],axis=1).ravel()-op.forward(pilot)
            bias=weights@expected-ell@delta;active=sd>1e-15
            coverage=(np.abs(bias)<=half+1e-12).astype(float)
            coverage[active]=norm.cdf((half[active]-bias[active])/sd[active])-norm.cdf((-half[active]-bias[active])/sd[active])
            true_target=float(ell@(pilot+delta));expected_center=true_target+bias
            positive=(expected_center>half).astype(float);negative=(expected_center < -half).astype(float)
            positive[active]=norm.cdf((expected_center[active]-half[active])/sd[active])
            negative[active]=norm.cdf((-expected_center[active]-half[active])/sd[active])
            for idx,name in enumerate(names):
                records.append({'truth':label,'method':name,'bias':float(bias[idx]),'analytic_coverage':float(coverage[idx]),
                                'true_target':true_target,'expected_center':float(expected_center[idx]),
                                'probability_positive_interval':float(positive[idx]),'probability_negative_interval':float(negative[idx]),
                                'correct_sign_probability':float((positive if true_target>0 else negative)[idx]) if true_target else None,
                                'wrong_sign_probability':float((negative if true_target>0 else positive)[idx]) if true_target else float(positive[idx]+negative[idx]),
                                'width_fraction_of_no_data':float(half[idx]/no_data),'density_bound_satisfied':bool(np.linalg.norm(delta)<=B+1e-12)})
        result['angles'].append({'angle_deg':angle_deg,'shift_pixels':shift,'pose_term_seconds':setup_seconds,'fits':fits,
                                'gaussian_prior_checks':gaussian_checks,'coverage':records})
    result['seconds']=time.perf_counter()-start
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--particles',type=int,default=128)
    p.add_argument('--box',type=int,default=24);p.add_argument('--frequency-radius',type=float,default=5.);p.add_argument('--support',type=float,default=.35)
    p.add_argument('--radius',type=float,default=2.);p.add_argument('--width',type=float,default=.07);p.add_argument('--offset',type=float,default=.08)
    p.add_argument('--target',choices=['contrast','center'],default='contrast')
    p.add_argument('--snr',type=float,default=.1);p.add_argument('--angles',default='.1,.5,2');p.add_argument('--shift',type=float,default=.01)
    p.add_argument('--maxiter',type=int,default=100);p.add_argument('--tolerance',type=float,default=.005)
    p.add_argument('--quadratic-pose',action='store_true');p.add_argument('--gaussian-baselines',action='store_true')
    p.add_argument('--seed',type=int,default=609315);p.add_argument('--signals',type=int,default=10)
    p.add_argument('--output',default='ambient-pose.json');args=p.parse_args()
    out=ROOT/'results/uncertainty/development'/args.output;result={'config':vars(args),'cases':[]}
    for dataset in args.datasets.split(','):
        result['cases'].append(run(args,dataset));out.write_text(json.dumps(result,indent=2)+'\n')
