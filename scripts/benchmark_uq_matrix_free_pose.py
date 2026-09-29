#!/usr/bin/env python3
"""Compare a matrix-free randomized pose upper audit with the saved dense one."""
import argparse
import hashlib
import json
import resource
import time
from pathlib import Path
import numpy as np
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_random_spectral import gaussian_power_upper
from fourier_splats.uq_continuous_moments import integrated_cubic_remainder
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--dataset',default='10028')
    parser.add_argument('--target',default='center')
    parser.add_argument('--source',default='continuous-quadrature-optimized')
    parser.add_argument('--angle',type=float,default=1.)
    parser.add_argument('--iterations',type=int,default=12)
    parser.add_argument('--probes',type=int,default=4)
    parser.add_argument('--failure-probability',type=float,default=1e-6)
    parser.add_argument('--seed',type=int,default=609591)
    parser.add_argument('--output',default='matrix-free-pose-probe')
    args=parser.parse_args();out=BASE/args.output;out.mkdir(exist_ok=True)
    path=out/f'{args.dataset}-{args.target}-{args.angle:g}.json'
    if path.exists():raise RuntimeError('Preserve previous results; use a new output directory')
    snapshot=source_snapshot(ROOT,Path(__file__))
    fits=json.loads((BASE/args.source/f'{args.dataset}.json').read_text());cfg=fits['config']
    fit=next(t['fit'] for t in fits['targets'] if t['target']==args.target and t['width_fraction_field']==.07)
    g=particle_geometry(ROOT,args.dataset,'inference_half0',radius=cfg['frequency_radius'],count=cfg['particles'],seed=cfg['seed'])
    weights_path=BASE/args.source/f'{args.dataset}-{args.target}-0.07-weights.npz'
    weights=np.load(weights_path);np.testing.assert_array_equal(weights['indices'],g['indices'])
    w=weights['weights'];noise=float(weights['noise_std']);angle=np.deg2rad(args.angle);shift=.01/24
    start=time.perf_counter()
    op=PolynomialPoseFieldOperator(g['k'],g['q'],g['ctf'],w,noise,angle,shift,order=32)
    scaling=op.establish_group_scaling(callback=lambda r:print(r,flush=True))
    scaling_seconds=time.perf_counter()-start
    tick=time.perf_counter()
    random_bound=gaussian_power_upper(op.spatial_gram,op.shape[0],args.failure_probability,args.probes,
                                      args.iterations,args.seed,scaling['quadrature_gram_trace'],
                                      callback=lambda r:print(r,flush=True))
    spectral_seconds=time.perf_counter()-tick
    bound=np.sqrt(scaling['sum_group_scales']*(random_bound['eigenvalue_upper']+scaling['integration_eigenvalue_pad']))
    bound=min(bound,scaling['triangle_field_bound'])
    remainder=integrated_cubic_remainder(g['k'],g['q'],g['ctf'],w,noise,angle,shift,2.,1.)
    bias=fit['bias']+3*bound+remainder['remainder_bias']
    half=bias_aware_half_width_stable(np.linalg.norm(w),bias,.05-args.failure_probability)
    dense=None
    if cfg['particles']==128:
        original=json.loads((BASE/'continuous-pose-large'/f'{args.dataset}.json').read_text())
        record=next(r for r in original['records'] if r['target']==args.target and r['angle_degrees']==args.angle)
        dense=record['audit']['joint_field']['spectral']
    result={'stage':'post-review development matrix-free audit; supplied simulation bounds; not experimental calibration',
            'config':vars(args),'source_snapshot':snapshot,'particles':len(g['k']),
            'weights_sha256':hashlib.sha256(weights_path.read_bytes()).hexdigest(),
            'matrix_free_shape':list(op.shape),'scaling_seconds':scaling_seconds,
            'spectral_seconds':spectral_seconds,'total_seconds':time.perf_counter()-start,
            'peak_resident_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'memory_units':'ru_maxrss bytes on macOS; KiB on Linux (this script currently benchmarked on macOS)',
            'spectral_probability_bound':random_bound,'group_scale_sum':scaling['sum_group_scales'],
            'integration_pad':scaling['integration_eigenvalue_pad'],'field_bound':float(bound),
            'dense_joint_field_bound':dense,'relative_field_bound':float(bound/dense) if dense else None,
            'density_bias':fit['bias'],'moment_remainder':remainder['remainder_bias'],
            'uncapped_relative_half_width':float(half/(2*fit['target_norm'])),
            'selected_relative_half_width':float(min(1,half/(2*fit['target_norm']))),
            'alpha_noise':.05-args.failure_probability,'alpha_numerical':args.failure_probability,
            'alpha_total':.05,'complete':True}
    path.write_text(json.dumps(result,indent=2)+'\n')
    print('DONE',result['total_seconds'],result['relative_field_bound'],flush=True)


if __name__=='__main__':main()
