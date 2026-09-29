#!/usr/bin/env python3
"""Explicit assumption violations for fixed development estimators.

These tests are diagnostics, not extensions of a theorem to an invalid model.
Gaussian/t coverage below is exact noise-marginal coverage for the stated
synthetic distribution. Experimental density has no available coverage label.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.stats import norm,t
from fourier_splats.physics import ctf
from fourier_splats.uq_data import particle_geometry,VoxelReference,VoxelObservationOperator,SupportedVoxelOperator
from fourier_splats.uq_nonlinear_stress import exact_pose_adjoint
from fourier_splats.uq_subspace import matrix_free_certificate
from fourier_splats.uncertainty import bias_aware_half_width
from audit_uq_grid_refinement import model,target
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development';MAPS={'10028':'2660','10049':'6487','10076':'8434'}


def run(args,dataset,name):
    begin=time.perf_counter();folder=BASE/f'pose-grid-{name}-0.5';saved=np.load(folder/f'{dataset}-coarse-weights.npz')
    source=json.loads((folder/f'{dataset}.json').read_text());cfg=source['config'];box=cfg['coarse_box'];n=cfg['particles']
    g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=n,seed=609315)
    np.testing.assert_array_equal(g['indices'],saved['indices']);ck=np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
    op,pilot,mask,noise=model(g,ck,box,source['noise_std_physical']);ell=saved['ell'];B=2.;angle=np.deg2rad(.5);shift=.01
    reference=VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map',box=box).volume.ravel()[mask];reference/=np.linalg.norm(reference)
    delta=reference-pilot;signal=op.forward(delta);m=2*op.op.q
    fixed=matrix_free_certificate(op,ell,B,rtol=.005)
    shared=json.loads((BASE/f'shared-density-{name}-0.5'/f'{dataset}.json').read_text())
    fits={'fixed_pose':(fixed.weights,fixed.half_width),
          'quadratic_block':(saved['weights'],source['coarse_half_width']),
          'quadratic_shared':(saved['weights'],shared['bounds']['minimum_valid_bound']['half_width'])}
    records=[]
    def record(axis,setting,mean,truth,noise_model='standard_Gaussian',sigma_multiplier=1.,rho=0.,df=None,signal_label='independent_map',**extra):
        for method,(w,half) in fits.items():
            bias=float(w@mean-truth);sd=np.linalg.norm(w);wb=w.reshape(n,m)
            if noise_model=='within_particle_equicorrelated':variance=(1-rho)*sd**2+rho*np.sum(wb.sum(axis=1)**2)
            elif noise_model=='global_weight_aligned_equicorrelated':variance=(1-rho)*sd**2+rho*np.sum(abs(w))**2
            else:variance=(sigma_multiplier*sd)**2
            actual_sd=np.sqrt(variance)
            if df is None:coverage=float(norm.cdf((half-bias)/actual_sd)-norm.cdf((-half-bias)/actual_sd))
            else:
                scale=actual_sd*np.sqrt((df-2)/df);coverage=float(t.cdf((half-bias)/scale,df)-t.cdf((-half-bias)/scale,df))
            records.append({'axis':axis,'setting':setting,'method':method,'signal':signal_label,'noise_model':noise_model,
                'analytic_coverage':coverage,'actual_bias':bias,'assumed_noise_sd':float(sd),'actual_noise_sd':float(actual_sd),
                'half_width':float(half),'df':df,**extra})
    truth=float(ell@delta)
    record('reference_control',0,signal,truth)
    # Each estimator's own nominal bias boundary is diagnostic, not a common truth.
    for owner,(w,_) in fits.items():
        h=op.adjoint(w)-ell;boundary=B*h/np.linalg.norm(h);mu=op.forward(boundary);tt=float(ell@boundary)
        record('nominal_density_boundary',1.,mu,tt,signal_label='boundary_for_'+owner)
        for multiplier in [1.5,2.,4.]:record('density_radius_violation',multiplier,mu*multiplier,tt*multiplier,signal_label='boundary_for_'+owner)
        for scale in [1.,1.5,2.,4.]:record('noise_scale',scale,mu,tt,sigma_multiplier=scale,signal_label='boundary_for_'+owner)
        for rho in [.1,.5,.9]:
            record('noise_correlation',rho,mu,tt,noise_model='within_particle_equicorrelated',rho=rho,signal_label='boundary_for_'+owner)
            record('adversarial_noise_correlation',rho,mu,tt,noise_model='global_weight_aligned_equicorrelated',rho=rho,signal_label='boundary_for_'+owner,
                   interpretation='artificial estimator-aligned global noise mode; not a typical experimental noise claim')
        for df in [3,5,10]:record('heavy_tails',df,mu,tt,noise_model='isotropic_multivariate_Student_t',df=df,signal_label='boundary_for_'+owner)
    # CTF errors affect both pilot and perturbation, while the estimator centers nominally.
    metadata=np.load(ROOT/'data'/dataset/'metadata.npz')['ctf'][g['indices']]
    for error in [-2000,-500,-100,0,100,500,2000]:
        params=metadata.copy();params[:,2:4]+=error;transfer=ctf(g['q'][0]/g['field_A'],params)
        actual=SupportedVoxelOperator(VoxelObservationOperator(g['k'],transfer,box,noise*box**1.5),mask)
        record('coherent_defocus_error_A',error,actual.forward(reference)-op.forward(pilot),truth)
    # Support-only violation: stay within the same full-grid L2 density radius.
    ellfull=target(box,np.ones(box**3,dtype=bool),.07,name);outside=~mask
    for owner,(w,_) in fits.items():
        hfull=op.op.adjoint(w)-ellfull;v=np.zeros(box**3);v[outside]=hfull[outside];v/=np.linalg.norm(v)
        remaining=np.sqrt(max(0,B**2-np.linalg.norm(delta)**2))
        for fraction in [.1,.5,1.]:
            amplitude=fraction*remaining;mu=signal+amplitude*op.op.forward(v);tt=truth+amplitude*float(ellfull@v)
            record('support_violation',fraction,mu,tt,signal_label='outside_direction_for_'+owner,
                   total_density_delta_norm=float(np.sqrt(np.linalg.norm(delta)**2+amplitude**2)),outside_density_norm=float(amplitude))
    # Two equally populated states have exactly the reference as their mean.
    direction=ell/np.linalg.norm(ell);projected=op.forward(direction).reshape(n,m)
    rng=np.random.default_rng(args.seed)
    for owner,(w,_) in fits.items():
        score=np.sum(w.reshape(n,m)*projected,axis=1);order=np.argsort(score);adversarial=np.full(n,-1.);adversarial[order[n//2:]]=1
        random=adversarial[rng.permutation(n)]
        for assignment,labels in [('balanced_random',random),('view_estimator_coupled',adversarial)]:
            assert labels.sum()==0
            for amplitude in [.1,.3,1.]:
                mean=signal+amplitude*(labels[:,None]*projected).ravel()
                record('two_state_heterogeneity',amplitude,mean,truth,signal_label=assignment+'_for_'+owner,
                       interpretation='estimand is the equal-population mean; states are reference plus/minus a target-directed density perturbation')
    # Exact rotations and shifts exceed (or obey) the prescribed pose ball.
    stressfolder=BASE/f'nonlinear-stress-{name}-0.5';stress=json.loads((stressfolder/f'{dataset}.json').read_text())
    winner=max(stress['candidates'],key=lambda r:r['absolute_bias']);adversary=np.load(stressfolder/f"{dataset}-sign{winner['sign']}-restart{winner['restart']}.npz")
    for multiplier in [0,1,2,4,8]:
        pose=adversary['pose']*multiplier;rho_true=pilot+adversary['density_delta'];tt=float(ell@adversary['density_delta'])
        # Use the adjoint identity separately for each estimator, retaining exact nonlinear physics.
        for method,(w,half) in fits.items():
            h=exact_pose_adjoint(op,g['q'],w,pose,angle,shift);bias=float(rho_true@h-pilot@op.adjoint(w)-tt);sd=np.linalg.norm(w)
            records.append({'axis':'joint_pose_radius','setting':multiplier,'method':method,'signal':'shared_estimator_nonlinear_adversary',
                'noise_model':'standard_Gaussian','analytic_coverage':float(norm.cdf((half-bias)/sd)-norm.cdf((-half-bias)/sd)),
                'actual_bias':bias,'assumed_noise_sd':float(sd),'actual_noise_sd':float(sd),'half_width':float(half),
                'maximum_normalized_pose_norm':float(np.linalg.norm(pose,axis=1).max())})
    out=BASE/args.output;out.mkdir(parents=True,exist_ok=True)
    result={'stage':'development assumption stresses; all settings retained; violations are outside the theorem',
            'dataset':dataset,'target':name,'config':vars(args),'records':records,'seconds':time.perf_counter()-begin,
            'fixed_pose_optimization_converged':fixed.converged,'source_weights_sha256':source['weights_sha256'],
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'interpretation':'Same fixed estimators; exact Gaussian or elliptical-t projected noise laws. Artificial adversaries must not be described as typical experimental cases.'}
    (out/f'{dataset}-{name}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(dataset,name,'records',len(records),'seconds',result['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--targets',default='center,contrast')
    p.add_argument('--seed',type=int,default=609360);p.add_argument('--output',default='assumption-stress');args=p.parse_args()
    for dataset in args.datasets.split(','):
        for name in args.targets.split(','):run(args,dataset,name)
