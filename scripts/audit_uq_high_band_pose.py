#!/usr/bin/env python3
"""Continuous pose post-audit of a completed higher-band fixed-pose fit.

Cheap quadrature chooses positive scales only. The high-order polynomial
operator, analytic integration pad and fresh spectral probes certify the bound.
No smaller design grid is substituted for the continuous density class.
"""
import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path
import numpy as np
import finufft
from fourier_splats.uq_data import particle_geometry,VoxelReference
from fourier_splats.uq_continuous import cell_forward,cell_target_coefficients
from fourier_splats.uq_continuous_pose import polynomial_kernel_error,pose_cell_forward
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_joint_bias import (pose_scale_record,scaled_pose_radius,
    sharp_cube_cubic_coefficients,residual_pose_cross_bound,joint_density_pose_bias,pair_pose_fourier_moments)
from fourier_splats.uq_cell_moments import cell_fourier_moments,check_cell_pose_pairings
from fourier_splats.uq_random_spectral import gaussian_power_upper
from fourier_splats.uq_intervals import bias_aware_half_width_stable,reference_interval_summary
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
MAPS={'10028':'2660','10049':'6487','10076':'8434'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--target',default='center')
    p.add_argument('--angle',type=float,default=1.);p.add_argument('--shift-A',type=float,default=.5)
    p.add_argument('--order',type=int,default=80);p.add_argument('--design-order',type=int,default=12)
    p.add_argument('--iterations',type=int,default=40);p.add_argument('--threads',type=int,default=2)
    p.add_argument('--certificate-seed',type=int,default=610271)
    p.add_argument('--alpha-total',type=float,default=.05);p.add_argument('--delta',type=float,default=1e-6)
    p.add_argument('--pilot-pairing',action='store_true')
    p.add_argument('--output',default='continuous-high-band-pose-1024-10A');args=p.parse_args()
    source=ROOT/args.fit;prior=json.loads(source.read_text())
    if not prior.get('complete') or prior.get('error'):raise ValueError('Completed fixed-pose fit required')
    if args.threads>1 and ('tmp/nufft-openmp' not in finufft.__file__ or 'torch' in sys.modules):
        raise RuntimeError('Threaded Mac run requires isolated CPU-only FINUFFT')
    rows=[r for r in prior['targets'] if r['target']==args.target]
    if len(rows)!=1:raise ValueError('Source must contain exactly one requested target/width')
    row=rows[0];fit=row['fit'];cfg=prior['config'];dataset=prior['dataset'];width=row['width_fraction_field']
    if not 0<args.delta<args.alpha_total<1:raise ValueError('Require 0 < delta < alpha_total < 1')
    if 'centers_fraction_field' in row:
        centers=row['centers_fraction_field'];signs=row['signs']
        if 'target_lock_sha256' in prior:
            from fourier_splats.uq_pilot_targets import read_target_lock,LOCK_SHA256
            lock=read_target_lock(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
            locked=next(d for d in lock['datasets'] if d['dataset']==dataset)
            feature=next(f for f in locked['features'] if f['name']==args.target)
            if prior['target_lock_sha256']!=LOCK_SHA256 or signs!=[1]:raise ValueError('Changed pilot-selected target')
            np.testing.assert_array_equal(centers,[feature['center_fraction_field']])
            np.testing.assert_allclose(width,locked['width_fraction_field'],rtol=1e-14)
            np.testing.assert_allclose([args.alpha_total,args.delta],[cfg['alpha_total'],cfg['delta']],rtol=1e-14)
    elif args.target in ['center','contrast']:
        centers=[[0,0,0]] if args.target=='center' else [[0,0,.08],[0,0,-.08]]
        signs=[1] if args.target=='center' else [1,-1]
    else:raise ValueError('Generic target requires explicit centers and signs')
    if 'pose_polynomial_bias' in fit or 'cubic_bias' in fit or 'fixed poses/CTFs' not in prior.get('stage',''):
        raise ValueError('Source bias must be a fixed-pose density residual only')
    out=BASE/args.output;out.mkdir(parents=True,exist_ok=True);path=out/f'{dataset}-{args.target}-{args.angle:g}.json'
    if path.exists():raise RuntimeError('Preserve previous outcome')
    result={'complete':False,'stage':'Conditional known-noise higher-band fixed-weight pose post-audit; not pose-optimized or experimentally calibrated',
        'config':vars(args),'dataset':dataset,'target':args.target,'width_fraction_field':width,
        'pose_set':'Per-particle joint ball: squared scaled rotation norm plus squared scaled translation norm <= 1',
        'source_fit_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'source_snapshot':source_snapshot(ROOT,Path(__file__),['scripts/audit_uq_grid_refinement.py'])}
    start=time.perf_counter()
    def save():path.write_text(json.dumps(result,indent=2)+'\n')
    def progress(message):print(message,flush=True)
    save()
    try:
        g=particle_geometry(ROOT,dataset,'inference_half0',radius=cfg['frequency_radius'],count=cfg['particles'],seed=cfg['seed'])
        wp=source.with_name(f'{dataset}-{args.target}-{width}-weights.npz');saved=np.load(wp)
        w=saved['weights'];noise=float(saved['noise_std']);np.testing.assert_array_equal(saved['indices'],g['indices'])
        result['source_weights_sha256']=hashlib.sha256(wp.read_bytes()).hexdigest()
        angle=np.deg2rad(args.angle);shift=args.shift_A/g['field_A'];B=2.;P=1.;delta=args.delta
        if prior.get('density_radius',2.)!=B or prior.get('pilot_norm_bound',1.)!=P:
            raise ValueError('This runner requires the declared B=2/P=1 class')
        result.update(centers_fraction_field=centers,signs=signs,target_lock_sha256=prior.get('target_lock_sha256'))
        op=PolynomialPoseFieldOperator(g['k'],g['q'],g['ctf'],w,noise,angle,shift,order=args.order,nthreads=args.threads)
        scaling=op.establish_quadrature_design_scaling(args.design_order,callback=progress)
        result['pose_scaling']=pose_scale_record(scaling['group_scales'],op.denominators)
        result['scaling_method']=scaling['scaling_method'];result['setup_seconds']=time.perf_counter()-start
        wr=w.reshape(op.n,2*op.nq);amp=np.hypot(wr[:,:op.nq],wr[:,op.nq:])
        masses=np.einsum('naj,nj->na',op.coefficient_bound_matrix(),amp)
        pad=float(polynomial_kernel_error(g['k'],args.order)*np.sum(masses*masses))
        save();spectral=gaussian_power_upper(op.spatial_gram,op.shape[0],delta,4,args.iterations,args.certificate_seed,callback=progress)
        L=scaled_pose_radius(scaling['group_scales']);field=L*np.sqrt(spectral['eigenvalue_upper']+pad)
        cubic=float(np.sum(sharp_cube_cubic_coefficients(g['k'],g['q'],g['ctf']/noise,angle,shift,B+P)*amp))
        cross=residual_pose_cross_bound(op,w,noise,centers,signs,width);vector=cross.pop('cross_vector')
        joint=joint_density_pose_bias(fit['bias']/B,field,cross['norm_upper'],L,B,P,cubic/(B+P))
        checkpoint=np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
        cellop,pilot,_,check_noise=model(g,checkpoint,24,noise=noise);pilot=cellop.expand(pilot)
        np.testing.assert_allclose(check_noise,noise);np.testing.assert_allclose(np.linalg.norm(pilot),P,rtol=1e-12)
        pilot_vector=None
        if args.pilot_pairing:
            pilot_moments=cell_fourier_moments(g['k'],pilot,24,nthreads=args.threads)
            pilot_vector=pair_pose_fourier_moments(op,pilot_moments)
            result['pilot_numerical_check']=check_cell_pose_pairings(op,pilot,24,pilot_vector,pilot_moments)
            if not result['pilot_numerical_check']['passed']:
                save();raise AssertionError('Pilot direct-sum check failed; preserve outcome')
            result['original_joint_bound']=joint
            joint=joint_density_pose_bias(fit['bias']/B,field,cross['norm_upper'],L,B,P,cubic/(B+P),float(np.linalg.norm(pilot_vector)))
        half=bias_aware_half_width_stable(fit['noise_sd'],joint['bias_upper'],args.alpha_total-delta);no_data=B*fit['target_norm']
        result.update(spectral_upper_certificate=spectral,pose_integration_pad=pad,cross=cross,bound=joint,
            density_radius=B,pilot_norm_bound=P,
            noise_sd=fit['noise_sd'],density_bias=fit['bias'],pose_polynomial_bias=(B+P)*field,cubic_bias=cubic,
            alpha_noise=args.alpha_total-delta,alpha_numerical=delta,alpha_total=args.alpha_total,
            half_width=half,selected_relative_half_width=min(1.,half/no_data),uses_no_data=bool(half>=no_data),
            particles=len(g['k']),fourier_pairs_per_particle=op.nq,target_width_A=width*g['field_A'],reference_checks=[])
        arrays={'residual_pose_cross':vector}
        if pilot_vector is not None:arrays['pilot_pose_cross']=pilot_vector
        np.savez(path.with_suffix('.npz'),**arrays)
        save()
        pilot_target=float(cell_target_coefficients(24,centers,signs,width)@pilot)
        result['pilot_target']=pilot_target
        nominal=cell_forward(g['k'],g['ctf'],pilot,24,noise)
        rho=VoxelReference.from_mrc(ROOT/f'data/uncertainty/references/emd_{MAPS[dataset]}.map',box=64).volume.ravel()
        rho/=np.linalg.norm(rho);truth=float(cell_target_coefficients(64,centers,signs,width)@rho)
        rng=np.random.default_rng(610281+int(dataset))
        for scenario in ['nominal','coherent_x','random_boundary']:
            poses=np.zeros((op.n,5))
            if scenario=='coherent_x':poses[:,0]=1.
            if scenario=='random_boundary':
                poses=rng.normal(size=poses.shape);poses/=np.linalg.norm(poses,axis=1)[:,None]
            signal=pose_cell_forward(g['k'],g['q'],g['ctf'],rho,64,noise,poses,angle,shift)
            expected=float(pilot_target+w@(signal-nominal))
            check=reference_interval_summary(truth,expected,fit['noise_sd'],pilot_target,half,no_data)
            result['reference_checks'].append({'scenario':scenario,'raw_expected_center':expected,**check});save()
        failure=any(r['analytic_coverage']<(1-args.alpha_total+delta)-1e-8 for r in result['reference_checks'])
        result.update(complete=True,reference_coverage_failure=failure,seconds=time.perf_counter()-start,
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
            minimum_reference_power=min(r['correct_sign_probability'] for r in result['reference_checks']))
        save()
        print('DONE',result['selected_relative_half_width'],result['minimum_reference_power'],flush=True)
    except Exception as exc:
        result.update(complete=False,error=repr(exc),seconds=time.perf_counter()-start);save();raise


if __name__=='__main__':main()
