#!/usr/bin/env python3
"""Actual large-design pose audit, separately labeled from useful reconstruction.

Uses distinct archived particle geometries, not replicated copies. The estimator
is a deterministic matched filter, not an optimized pose-aware reconstruction.
Its possible vacuity is retained. No new particle pixels are needed for this
conditional numerical scalability experiment.
"""
import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path
import numpy as np
from fourier_splats.physics import ctf
from fourier_splats.uq_data import half_plane_frequencies
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_pose_optimization import cubic_penalty_coefficients
from fourier_splats.uq_random_spectral import gaussian_power_upper
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    p=argparse.ArgumentParser();p.add_argument('--dataset',default='10028')
    p.add_argument('--particles',type=int,default=10000);p.add_argument('--radius',type=float,default=12)
    p.add_argument('--order',type=int,default=80);p.add_argument('--iterations',type=int,default=24)
    p.add_argument('--angle',type=float,default=1.);p.add_argument('--shift-A',type=float,default=.5)
    p.add_argument('--output',default='matrix-free-pose-10000-radius12');args=p.parse_args()
    out=BASE/args.output;out.mkdir(exist_ok=True);path=out/f'{args.dataset}.json'
    if path.exists():raise RuntimeError('Preserve existing run; choose another output')
    snapshot=source_snapshot(ROOT,Path(__file__));start=time.perf_counter()
    paths=[ROOT/f'data/{args.dataset}',ROOT/f'data/uncertainty/confirmation/prediction-v1/{args.dataset}']
    metadata=[np.load(path/'metadata.npz') for path in paths]
    indices=np.concatenate([m['indices'] for m in metadata])
    if len(np.unique(indices))!=len(indices):raise AssertionError('Expected distinct archive particle indices')
    if args.particles>len(indices):raise ValueError('Insufficient distinct available geometry')
    selection=np.random.default_rng(609641).permutation(len(indices))[:args.particles]
    rotations=np.concatenate([m['rotations'] for m in metadata])[selection]
    parameters=np.concatenate([m['ctf'] for m in metadata])[selection]
    manifest=json.loads((paths[0]/'manifest.json').read_text())
    field=manifest['raw_box']*manifest['raw_pixel_size_A']
    q0=half_plane_frequencies(args.radius);q=np.broadcast_to(q0,(args.particles,len(q0),2))
    plane=np.pad(q0,((0,0),(0,1)));k=np.einsum('qi,nij->nqj',plane,rotations)
    transfer=ctf(q0/field,parameters).astype(float)
    noise=float(np.load(BASE/'continuous-quadrature-optimized'/f'{args.dataset}-center-0.07-weights.npz')['noise_std'])
    gram=QuadratureObservationGram(k,transfer,noise,order=args.order,preconditioner_rank=0)
    a,norm2=gram.target([[0,0,0]],[1],.07)
    w=norm2*a/(a@a)
    gw=gram.matvec(w);error=gram.quadrature_error(w)
    residual=np.sqrt(max(0,norm2-2*w@a+w@gw+error['squared_field_norm']))
    del gw,a,gram
    angle=np.deg2rad(args.angle);shift=args.shift_A/field;B=2.;P=1.;delta=1e-6
    op=PolynomialPoseFieldOperator(k,q,transfer,w,noise,angle,shift,order=args.order)
    scaling=op.establish_coefficient_scaling()
    result={'stage':'Post-review conditional numerical scale experiment; distinct experimental geometries and deterministic matched-filter weights; not an optimized reconstruction or empirical calibration',
            'config':vars(args),'source_snapshot':snapshot,'complete':False,
            'particles':args.particles,'independent_fourier_pairs_per_particle':len(q0),
            'archive_indices_sha256':hashlib.sha256(indices[selection].tobytes()).hexdigest(),
            'field_A':field,'noise_std_supplied':noise,'target_width_fraction_field':.07,
            'density_radius':B,'pilot_norm_bound':P,'rotation_radius_degrees':args.angle,
            'translation_radius_A':args.shift_A,'setup_seconds':time.perf_counter()-start,
            'scaling_method':scaling['scaling_method'],'matrix_free_shape':list(op.shape),
            'power_history':[]}
    np.savez(out/f'{args.dataset}-weights.npz',weights=w,archive_indices=indices[selection],noise_std=noise)
    path.write_text(json.dumps(result,indent=2)+'\n');print('SETUP',result['setup_seconds'],op.shape,flush=True)
    def progress(row):
        result['power_history'].append(row);result['elapsed_seconds']=time.perf_counter()-start
        path.write_text(json.dumps(result,indent=2)+'\n');print(row,flush=True)
    try:
        spectral=gaussian_power_upper(op.spatial_gram,op.shape[0],delta,4,args.iterations,609651,
                                      scaling['quadrature_gram_trace_upper'],callback=progress)
        field_bound=np.sqrt(scaling['sum_group_scales']*(spectral['eigenvalue_upper']+scaling['integration_eigenvalue_pad']))
        field_bound=min(field_bound,scaling['triangle_field_bound'])
        amplitudes=np.hypot(w.reshape(args.particles,-1)[:,:len(q0)],w.reshape(args.particles,-1)[:,len(q0):])
        remainder=float(np.sum(cubic_penalty_coefficients(k,q,transfer/noise,angle,shift,B+P)*amplitudes))
        bias=B*residual+(B+P)*field_bound+remainder
        half=bias_aware_half_width_stable(np.linalg.norm(w),bias,.05-delta)
        result.update(complete=True,total_seconds=time.perf_counter()-start,
                      peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
                      density_bias=float(B*residual),pose_polynomial_bias=float((B+P)*field_bound),
                      cubic_bias=remainder,uncapped_relative_half_width=float(half/(B*np.sqrt(norm2))),
                      selected_relative_half_width=float(min(1,half/(B*np.sqrt(norm2)))),
                      spectral_upper_certificate=spectral,alpha_noise=.05-delta,alpha_numerical=delta,alpha_total=.05)
    except Exception as exc:
        result.update(error=repr(exc),elapsed_seconds=time.perf_counter()-start)
        path.write_text(json.dumps(result,indent=2)+'\n');raise
    path.write_text(json.dumps(result,indent=2)+'\n');print('DONE',result['total_seconds'],result['selected_relative_half_width'],flush=True)


if __name__=='__main__':main()
