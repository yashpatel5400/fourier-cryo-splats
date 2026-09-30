#!/usr/bin/env python3
"""Fit with a training covariance proxy; audit noise on independent exposures."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from fourier_splats.physics import fft_center
from fourier_splats.uq_data import particle_geometry,particle_observations
from fourier_splats.uq_continuous import continuous_certificate,cell_forward,cell_target_coefficients
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_continuous_pose import continuous_pose_audit
from fourier_splats.uq_continuous_moments import integrated_cubic_remainder
from fourier_splats.uq_noise_design import CenteredNoiseMetric,NoiseMetricGram,shrunk_second_moment
from fourier_splats.uq_noise_calibration import estimator_variance_upper,uncenter_fourier_weights
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def realify(x):return np.concatenate([x.real,x.imag],axis=1)


def run(args,dataset,snapshot):
    source=BASE/'experimental-noise-grouped'/f'{dataset}.json';old=json.loads(source.read_text())
    if not old['complete']:raise RuntimeError('Completed source cohort required')
    out=BASE/args.output;out.mkdir(exist_ok=True);path=out/f'{dataset}.json'
    if path.exists():raise RuntimeError('Preserve earlier results')
    start=time.perf_counter();g=particle_geometry(ROOT,dataset,'inference_half0',radius=5)
    lookup={int(index):i for i,index in enumerate(g['indices'])};chosen=np.array([lookup[index] for index in old['inference_indices']])
    for key in ['indices','groups','q','k','ctf','translations']:g[key]=g[key][chosen]
    checkpoint=np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
    op,pilot,_,_=model(g,checkpoint,24);pilot=op.expand(pilot)
    amplitude=old['pilot_amplitude_raw_units'];images=np.load(ROOT/'data/uncertainty/confirmation/prediction-v1'/dataset/'images.npy',mmap_mode='r')
    representatives=np.array(old['representative_indices']);split=old['calibration_groups'];q=g['q'][0].astype(int);center=images.shape[-1]//2
    fourier=fft_center(np.asarray(images[representatives],dtype=float));samples=realify(fourier[:,q[:,1]+center,q[:,0]+center])/amplitude
    training,calibration=samples[:split],samples[split:]
    covariance=shrunk_second_moment(training,args.shrinkage)
    phases=g['manifest']['data_sign']*np.exp(-2j*np.pi*np.einsum('nqi,ni->nq',g['q'],g['translations']))
    metric=CenteredNoiseMetric(covariance,phases)
    gram=QuadratureObservationGram(g['k'],g['ctf'],1.,order=40,preconditioner_rank=args.preconditioner_rank)
    transformed=NoiseMetricGram(gram,metric)
    y=realify(particle_observations(ROOT,dataset,g)).ravel()/amplitude
    nominal=cell_forward(g['k'],g['ctf'],pilot,24,1.);B=2.;P=float(np.linalg.norm(pilot));alpha=.045
    result={'stage':'Exploratory covariance-guided design with independently audited noise; common-noise, pose and density assumptions unverified',
            'dataset':dataset,'config':vars(args),'source_snapshot':snapshot,'complete':False,
            'source_result_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'inference_indices':g['indices'].tolist(),
            'training_groups':len(training),'calibration_groups':len(calibration),
            'calibration_representative_indices':representatives[split:].tolist(),
            'training_metric_eigenvalue_range':[float(metric.values[0]),float(metric.values[-1])],
            'metric_scope':'Shrunk uncentered empirical second moment is only an optimization guide, not a covariance confidence bound.',
            'noise_alpha':alpha,'calibration_alpha':.005,'unverified_conditions':old['unverified_conditions'],
            'selection_scope':'Exploratory post-review development, not frozen confirmation or simultaneous coverage.',
            'records':[]}
    def save():path.write_text(json.dumps(result,indent=2)+'\n')
    save()
    try:
        for width in sorted(set(row['target_width_fraction'] for row in old['records']),reverse=True):
            def progress(row):print('FIT',dataset,width,row['iteration'],row['relative_gap'],row['cg_iterations'],flush=True)
            fit=continuous_certificate(transformed,[[0,0,0]],[1],width,B,alpha=alpha,rtol=.005,maxiter=100,callback=progress)
            x=fit.pop('weights');w=metric.apply(x);raw=uncenter_fourier_weights(w.reshape(len(g['indices']),-1),phases)
            bound=estimator_variance_upper(calibration,raw,.005);sd=bound['noise_sd_upper']
            # Preserve transformed and physical weights and the learned metric.
            np.savez(out/f'{dataset}-width{width:.8f}-weights.npz',weights=w,design_coordinates=x,covariance_proxy=covariance,indices=g['indices'],noise_std=1.)
            pilot_target=float(cell_target_coefficients(24,[[0,0,0]],[1],width)@pilot)
            rawcenter=float(pilot_target+w@(y-nominal));no_data=B*fit['target_norm']
            oldrows=[r for r in old['records'] if r['target_width_fraction']==width]
            for oldrow in oldrows:
                degrees=oldrow['rotation_radius_degrees'];audit=None;moment=None
                if degrees:
                    audit=continuous_pose_audit(g['k'],g['q'],g['ctf'],w,1.,np.deg2rad(degrees),.5/g['field_A'],B,P,fit['bias']/B,order=32)
                    moment=integrated_cubic_remainder(g['k'],g['q'],g['ctf'],w,1.,np.deg2rad(degrees),.5/g['field_A'],B,P)
                    bias=audit['density_bias']+audit['pose_polynomial_bias']+min(audit['remainder_bias'],moment['remainder_bias'])
                else:bias=fit['bias']
                half=bias_aware_half_width_stable(sd,bias,alpha);fallback=bool(half>=no_data);chosencenter=pilot_target if fallback else rawcenter
                row={'sigma_A':oldrow['sigma_A'],'rotation_radius_degrees':degrees,'translation_radius_A':oldrow['translation_radius_A'],
                     'design_proxy_fit':fit,'directional_calibration':bound,'density_bias':fit['bias'],
                     'total_bias_upper':float(bias),'half_width_before_fallback':half,'relative_half_width':min(half,no_data)/no_data,
                     'interval_center':chosencenter,'interval_half_width':min(half,no_data),'uses_no_data':fallback,
                     'excludes_zero':bool(abs(chosencenter)>min(half,no_data)),
                     'approximate_reference_target':oldrow['approximate_reference_target'],
                     'approximate_reference_inside_interval':bool(abs(oldrow['approximate_reference_target']-chosencenter)<=min(half,no_data))}
                result['records'].append(row);save();print('CASE',dataset,width,degrees,row['relative_half_width'],row['excludes_zero'],flush=True)
        result.update(complete=True,seconds=time.perf_counter()-start);save();print('DONE',dataset,result['seconds'],flush=True)
    except Exception as exc:
        result.update(error=repr(exc),seconds=time.perf_counter()-start);save();raise


def main():
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--output',default='noise-metric-design')
    p.add_argument('--shrinkage',type=float,default=.2);p.add_argument('--preconditioner-rank',type=int,default=256)
    args=p.parse_args();snapshot=source_snapshot(ROOT,Path(__file__),['scripts/audit_uq_grid_refinement.py'])
    for dataset in args.datasets.split(','):run(args,dataset,snapshot)


if __name__=='__main__':main()
