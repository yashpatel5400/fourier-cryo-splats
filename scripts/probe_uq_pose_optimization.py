#!/usr/bin/env python3
"""Bounded development probe of continuous pose-aware estimator weights."""
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

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def main():
    p=argparse.ArgumentParser();p.add_argument('--dataset',default='10028')
    p.add_argument('--target',default='center');p.add_argument('--angle',type=float,default=1.)
    p.add_argument('--maxiter',type=int,default=20);p.add_argument('--power-iterations',type=int,default=24)
    p.add_argument('--output',default='pose-aware-optimization-probe');args=p.parse_args()
    out=BASE/args.output;out.mkdir(exist_ok=True);path=out/f'{args.dataset}-{args.target}-{args.angle:g}.json'
    if path.exists():raise RuntimeError('Preserve existing results')
    snapshot=source_snapshot(ROOT,Path(__file__))
    base=json.loads((BASE/'continuous-quadrature-optimized'/f'{args.dataset}.json').read_text())
    cfg=base['config'];g=particle_geometry(ROOT,args.dataset,'inference_half0',radius=5,count=128,seed=cfg['seed'])
    weight_path=BASE/'continuous-quadrature-optimized'/f'{args.dataset}-{args.target}-0.07-weights.npz'
    saved=np.load(weight_path);np.testing.assert_array_equal(saved['indices'],g['indices'])
    centers=[[0,0,0]] if args.target=='center' else [[0,0,.08],[0,0,-.08]]
    signs=[1] if args.target=='center' else [1,-1]
    start=time.perf_counter()
    result={'stage':'Exploratory post-review continuous pose-aware optimization; independent conic validation still required',
            'config':vars(args),'source_snapshot':snapshot,
            'initial_weights_sha256':hashlib.sha256(weight_path.read_bytes()).hexdigest(),'complete':False}
    path.write_text(json.dumps(result,indent=2)+'\n')
    gram=QuadratureObservationGram(g['k'],g['ctf'],float(saved['noise_std']),order=40,preconditioner_rank=1024)
    try:
        fit=pose_aware_certificate(gram,g['q'],centers,signs,.07,2.,1.,saved['weights'],
                                  np.deg2rad(args.angle),.01/24,maxiter=args.maxiter,
                                  power_iterations=args.power_iterations,callback=lambda r:print(r,flush=True))
        w=fit.pop('weights');np.savez(out/f'{args.dataset}-{args.target}-{args.angle:g}-weights.npz',weights=w,
                                     indices=g['indices'],noise_std=float(saved['noise_std']))
        result.update(fit=fit,complete=True,seconds=time.perf_counter()-start)
        result['relative_half_width']=min(1.,fit['half_width']/(2*fit['target_norm']))
    except Exception as exc:
        result.update(error=repr(exc),seconds=time.perf_counter()-start)
        path.write_text(json.dumps(result,indent=2)+'\n')
        raise
    path.write_text(json.dumps(result,indent=2)+'\n')
    print('DONE',result['relative_half_width'],fit['relative_sum_gap'],flush=True)


if __name__=='__main__':main()
