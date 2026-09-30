#!/usr/bin/env python3
"""Numerical conformance and timing of transform-plan reuse on real geometry."""
import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path
import numpy as np
import finufft
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_planned_transforms import PlannedQuadratureObservationGram,PlannedPolynomialPoseFieldOperator
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    p=argparse.ArgumentParser();p.add_argument('--dataset',default='10049');p.add_argument('--particles',type=int,default=128)
    p.add_argument('--radius',type=float,default=5);p.add_argument('--order',type=int,default=32)
    p.add_argument('--repeats',type=int,default=3);p.add_argument('--threads',type=int,default=2)
    p.add_argument('--output',default='transform-plan-benchmark-128');args=p.parse_args()
    if args.threads>1 and sys.platform=='darwin' and ('tmp/nufft-openmp' not in finufft.__file__ or 'torch' in sys.modules):
        raise RuntimeError('Threaded Mac run requires isolated CPU-only FINUFFT')
    out=BASE/args.output;out.mkdir(exist_ok=True);path=out/f'{args.dataset}.json'
    if path.exists():raise RuntimeError('Preserve previous benchmark')
    result={'complete':False,'config':vars(args),'source_snapshot':source_snapshot(ROOT,Path(__file__)),
            'scope':'Implementation/timing control using real acquisition coordinates and seeded random coefficients; not an inference result. Concurrent Mac load, not isolated timings.',
            'records':[]};start=time.perf_counter()
    def save():path.write_text(json.dumps(result,indent=2)+'\n')
    save()
    try:
        g=particle_geometry(ROOT,args.dataset,'inference_half0',radius=args.radius,count=args.particles,seed=609315)
        result['indices_sha256']=hashlib.sha256(g['indices'].tobytes()).hexdigest()
        rng=np.random.default_rng(610323);w=rng.normal(size=2*g['ctf'].size);w/=np.linalg.norm(w)
        for kind,classes in [('gram',[QuadratureObservationGram,PlannedQuadratureObservationGram]),
                             ('pose',[PolynomialPoseFieldOperator,PlannedPolynomialPoseFieldOperator])]:
            expected=None;reference_times=None
            for cls in classes:
                if kind=='gram':
                    op=cls(g['k'],g['ctf'],.02,order=args.order,preconditioner_rank=0);op.nthreads=args.threads
                    vector=w;function=op.matvec
                else:
                    op=cls(g['k'],g['q'],g['ctf'],w,.02,np.deg2rad(1.),.5/g['field_A'],order=args.order,nthreads=args.threads)
                    vector=np.random.default_rng(610324).normal(size=op.shape[0]);vector/=np.linalg.norm(vector);function=op.spatial_gram
                times=[]
                for trial in range(args.repeats):
                    begin=time.perf_counter();value=function(vector);times.append(time.perf_counter()-begin)
                    print(kind,cls.__name__,trial,times[-1],flush=True)
                    if expected is None:expected=value.copy()
                    np.testing.assert_allclose(value,expected,rtol=1e-9,atol=1e-10)
                difference=float(np.linalg.norm(value-expected)/max(np.linalg.norm(expected),1e-300))
                if reference_times is None:reference_times=times
                row={'operator':kind,'implementation':cls.__name__,'seconds':times,
                     'relative_output_difference':difference,
                     'warm_speed_ratio':float(np.mean(reference_times[1:])/np.mean(times[1:])) if len(times)>1 else None}
                result['records'].append(row);save();del op,function
        result.update(complete=True,seconds=time.perf_counter()-start,
                peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024));save()
    except Exception as exc:
        result.update(complete=False,error=repr(exc),seconds=time.perf_counter()-start);save();raise


if __name__=='__main__':main()
