#!/usr/bin/env python3
"""Post-review development grid of pose-aware fits with realistic shifts.

All cases and failures are retained. This is exploratory algorithm development,
not the frozen confirmation or an experimental nuisance-calibration claim.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_pose_optimization import pose_aware_certificate
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def run(args,dataset,snapshot):
    base_path=BASE/'continuous-quadrature-optimized'/f'{dataset}.json'
    base=json.loads(base_path.read_text());cfg=base['config']
    g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=128,seed=cfg['seed'])
    out=BASE/args.output;out.mkdir(exist_ok=True)
    saved=np.load(BASE/'continuous-quadrature-optimized'/f'{dataset}-center-0.07-weights.npz')
    noise=float(saved['noise_std'])
    gram=QuadratureObservationGram(g['k'],g['ctf'],noise,order=40,preconditioner_rank=1024)
    for target in args.targets.split(','):
        centers=[[0,0,0]] if target=='center' else [[0,0,.08],[0,0,-.08]]
        signs=[1] if target=='center' else [1,-1]
        wp=BASE/'continuous-quadrature-optimized'/f'{dataset}-{target}-0.07-weights.npz'
        saved=np.load(wp);np.testing.assert_array_equal(saved['indices'],g['indices'])
        for degrees in map(float,args.angles.split(',')):
            name=f'{dataset}-{target}-{degrees:g}';path=out/f'{name}.json'
            if path.exists():raise RuntimeError('Preserve prior cases; choose another output directory')
            start=time.perf_counter();seed=609701+int(dataset)*100+int(degrees*10)+(10000000 if target!='center' else 0)
            result={'stage':'Post-review pose-aware development; simulation noise and supplied nuisance bounds',
                    'dataset':dataset,'target':target,'rotation_radius_degrees':degrees,
                    'translation_radius_A':args.shift_A,'field_A':g['field_A'],'density_radius':2.,'pilot_norm_bound':1.,
                    'config':vars(args),'source_snapshot':snapshot,'source_geometry_config':cfg,
                    'initial_weights_sha256':hashlib.sha256(wp.read_bytes()).hexdigest(),'complete':False}
            path.write_text(json.dumps(result,indent=2)+'\n')
            def checkpoint(weights,row):
                np.savez(out/f'{name}-checkpoint.npz',weights=weights,indices=g['indices'],noise_std=noise)
                (out/f'{name}-checkpoint.json').write_text(json.dumps(row,indent=2)+'\n')
            def progress(row):print(name,row,flush=True)
            try:
                fit=pose_aware_certificate(gram,g['q'],centers,signs,.07,2.,1.,saved['weights'],
                                           np.deg2rad(degrees),args.shift_A/g['field_A'],maxiter=args.maxiter,
                                           power_iterations=args.power_iterations,certificate_seed=seed,
                                           callback=progress,checkpoint_callback=checkpoint)
                w=fit.pop('weights');np.savez(out/f'{name}-weights.npz',weights=w,indices=g['indices'],noise_std=noise)
                fallback=bool(fit['half_width']>=2*fit['target_norm'])
                result.update(fit=fit,complete=True,seconds=time.perf_counter()-start,
                              uses_no_data=fallback,uncapped_relative_half_width=fit['half_width']/(2*fit['target_norm']),
                              selected_relative_half_width=min(1.,fit['half_width']/(2*fit['target_norm'])))
            except Exception as exc:
                result.update(error=repr(exc),seconds=time.perf_counter()-start)
                path.write_text(json.dumps(result,indent=2)+'\n')
                print('FAILED_CASE',name,repr(exc),flush=True)
                continue
            path.write_text(json.dumps(result,indent=2)+'\n')
            print('COMPLETE_CASE',name,result['selected_relative_half_width'],fit['relative_sum_gap'],flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076')
    p.add_argument('--targets',default='center,contrast');p.add_argument('--angles',default='.5,1,2')
    p.add_argument('--shift-A',type=float,default=.5);p.add_argument('--maxiter',type=int,default=100)
    p.add_argument('--power-iterations',type=int,default=40);p.add_argument('--output',default='pose-aware-optimized-shift05')
    args=p.parse_args();snapshot=source_snapshot(ROOT,Path(__file__))
    if any(t not in ['center','contrast'] for t in args.targets.split(',')):raise ValueError('Unknown target')
    for dataset in args.datasets.split(','):run(args,dataset,snapshot)


if __name__=='__main__':main()
