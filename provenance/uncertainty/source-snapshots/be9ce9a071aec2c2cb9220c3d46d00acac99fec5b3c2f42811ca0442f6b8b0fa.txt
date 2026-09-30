#!/usr/bin/env python3
"""Target-specific, independently calibrated noise audit of existing estimators.

The first fresh-exposure half already selected the old fixed weights. Only the
second half calibrates their scalar noise variance here. Pull weights back to
raw Fourier coordinates before applying the common-covariance trace lemma.
No weights, targets, density classes, or original records are overwritten.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from fourier_splats.physics import fft_center
from fourier_splats.uq_data import particle_geometry,particle_observations
from fourier_splats.uq_continuous import cell_forward,cell_target_coefficients
from fourier_splats.uq_continuous_pose import continuous_pose_audit
from fourier_splats.uq_continuous_moments import integrated_cubic_remainder
from fourier_splats.uq_noise_calibration import common_covariance_trace_upper,estimator_variance_upper,uncenter_fourier_weights
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def realify(x):return np.concatenate([x.real,x.imag],axis=1)


def run(args,dataset,snapshot):
    source=BASE/'experimental-noise-grouped'/f'{dataset}.json';original=json.loads(source.read_text())
    if not original.get('complete'):raise ValueError('Completed source study required')
    out=BASE/args.output;out.mkdir(exist_ok=True);path=out/f'{dataset}.json'
    if path.exists():raise RuntimeError('Preserve all previous outcomes')
    start=time.perf_counter();g=particle_geometry(ROOT,dataset,'inference_half0',radius=5)
    lookup={int(index):i for i,index in enumerate(g['indices'])}
    chosen=np.array([lookup[index] for index in original['inference_indices']])
    for key in ['indices','groups','q','k','ctf','translations']:g[key]=g[key][chosen]
    checkpoint=np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
    op,pilot,_,_=model(g,checkpoint,24);pilot=op.expand(pilot)
    amplitude=original['pilot_amplitude_raw_units'];oldnoise=original['noise_sd_envelope']
    images=np.load(ROOT/'data/uncertainty/confirmation/prediction-v1'/dataset/'images.npy',mmap_mode='r')
    representatives=np.array(original['representative_indices']);split=original['calibration_groups']
    q=g['q'][0].astype(int);center=images.shape[-1]//2
    transformed=fft_center(np.asarray(images[representatives],dtype=float))
    samples=realify(transformed[:,q[:,1]+center,q[:,0]+center])/amplitude
    oldcal=common_covariance_trace_upper(samples[:split],.005)
    np.testing.assert_allclose(oldcal['covariance_trace_upper'],oldnoise**2,rtol=1e-12)
    # The original first calibration half may choose weights. The independently
    # held-out second half is the sole source of the NEW feature variance bound.
    calibration=samples[split:]
    y=realify(particle_observations(ROOT,dataset,g)).ravel()/(amplitude*oldnoise)
    nominal=cell_forward(g['k'],g['ctf'],pilot,24,oldnoise)
    phases=g['manifest']['data_sign']*np.exp(-2j*np.pi*np.einsum('nqi,ni->nq',g['q'],g['translations']))
    B=original['density_radius_supplied'];P=float(np.linalg.norm(pilot));alpha=.045
    result={'stage':'Exploratory independent directional noise calibration of existing fixed estimators; no experimental density coverage claim',
            'dataset':dataset,'config':vars(args),'source_snapshot':snapshot,'complete':False,
            'source_result':str(source.relative_to(ROOT)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'weight_selection_calibration_groups':split,'new_calibration_groups':len(calibration),
            'calibration_representative_indices':representatives[split:].tolist(),
            'inference_particles':len(g['indices']),'inference_indices':g['indices'].tolist(),
            'noise_alpha':alpha,'calibration_alpha':.005,'guarantee_scope':'Pointwise per fixed estimator, not simultaneous across targets or sensitivities',
            'unverified_conditions':original['unverified_conditions'],
            'raw_frame_note':'Calibration uses uncentered image Fourier coordinates; inference weights include the transpose of each supplied centering phase and the original whitening scale.',
            'post_hoc_note':'The second exposure half was previously used only for aggregate energy diagnostics; this development analysis is not a new frozen confirmation.',
            'records':[]}
    path.write_text(json.dumps(result,indent=2)+'\n')
    try:
        for width in sorted(set(row['target_width_fraction'] for row in original['records']),reverse=True):
            rows=[row for row in original['records'] if row['target_width_fraction']==width]
            weight_path=source.parent/f'{dataset}-width{width:.8f}-weights.npz';saved=np.load(weight_path)
            np.testing.assert_array_equal(saved['indices'],g['indices']);w=saved['weights']
            rawweights=uncenter_fourier_weights(w.reshape(len(g['indices']),-1),phases,oldnoise)
            bound=estimator_variance_upper(calibration,rawweights,.005)
            same_data_trace=common_covariance_trace_upper(calibration,.005)['covariance_trace_upper']*np.sum(rawweights**2)
            assert bound['estimator_variance_upper']<=same_data_trace*(1+1e-12)
            sd=bound['noise_sd_upper'];oldsd=float(np.linalg.norm(w))
            pilot_target=float(cell_target_coefficients(24,[[0,0,0]],[1],width)@pilot)
            rawcenter=float(pilot_target+w@(y-nominal));fit=rows[0]['fit'];no_data=B*fit['target_norm']
            for oldrow in rows:
                degrees=oldrow['rotation_radius_degrees'];audit=moment=None
                if degrees:
                    audit=continuous_pose_audit(g['k'],g['q'],g['ctf'],w,oldnoise,np.deg2rad(degrees),.5/g['field_A'],B,P,fit['bias']/B,order=32)
                    moment=integrated_cubic_remainder(g['k'],g['q'],g['ctf'],w,oldnoise,np.deg2rad(degrees),.5/g['field_A'],B,P)
                    bias=audit['density_bias']+audit['pose_polynomial_bias']+min(audit['remainder_bias'],moment['remainder_bias'])
                else:bias=fit['bias']
                oldhalf=bias_aware_half_width_stable(oldsd,bias,alpha)
                np.testing.assert_allclose(oldhalf,oldrow['half_width_before_fallback'],rtol=1e-9)
                half=bias_aware_half_width_stable(sd,bias,alpha);fallback=bool(half>=no_data)
                chosenhalf=min(half,no_data);chosencenter=pilot_target if fallback else rawcenter
                row={'sigma_A':oldrow['sigma_A'],'rotation_radius_degrees':degrees,
                     'translation_radius_A':oldrow['translation_radius_A'],'bias_upper_unchanged':float(bias),
                     'directional_calibration':bound,'old_noise_sd_upper':oldsd,
                     'noise_sd_ratio':sd/oldsd,'same_calibration_trace_sd_upper':float(np.sqrt(same_data_trace)),
                     'half_width_before_fallback':half,'old_half_width_before_fallback':oldhalf,
                     'relative_half_width':chosenhalf/no_data,'old_relative_half_width':oldrow['relative_half_width'],
                     'interval_center':chosencenter,'interval_half_width':chosenhalf,'uses_no_data':fallback,
                     'excludes_zero':bool(abs(chosencenter)>chosenhalf),
                     'approximate_reference_target':oldrow['approximate_reference_target'],
                     'approximate_reference_inside_interval':bool(abs(oldrow['approximate_reference_target']-chosencenter)<=chosenhalf),
                     'weights_sha256':hashlib.sha256(weight_path.read_bytes()).hexdigest()}
                result['records'].append(row);path.write_text(json.dumps(result,indent=2)+'\n')
                print('CASE',dataset,row['sigma_A'],degrees,'sd_ratio',row['noise_sd_ratio'],'relative_width',row['relative_half_width'],'excludes_zero',row['excludes_zero'],flush=True)
        result.update(complete=True,seconds=time.perf_counter()-start)
    except Exception as exc:
        result.update(error=repr(exc),seconds=time.perf_counter()-start);path.write_text(json.dumps(result,indent=2)+'\n');raise
    path.write_text(json.dumps(result,indent=2)+'\n');print('DONE',dataset,result['seconds'],flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--output',default='directional-noise-audit')
    args=p.parse_args();snapshot=source_snapshot(ROOT,Path(__file__),['scripts/audit_uq_grid_refinement.py'])
    for dataset in args.datasets.split(','):run(args,dataset,snapshot)


if __name__=='__main__':main()
