#!/usr/bin/env python3
"""Measure reusable FINUFFT plans before changing the pose implementation."""
import json
import time
from pathlib import Path
import numpy as np
import finufft
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def main():
    path = BASE/'nufft-plan-reuse-profile.json'
    if path.exists(): raise RuntimeError('Preserve the previous timing experiment')
    source = source_snapshot(ROOT, Path(__file__))
    g = particle_geometry(ROOT, '10028', 'inference_half0', radius=5, count=128, seed=609315)
    saved = np.load(BASE/'continuous-quadrature-optimized/10028-center-0.07-weights.npz')
    np.testing.assert_array_equal(saved['indices'], g['indices'])
    op = PolynomialPoseFieldOperator(g['k'], g['q'], g['ctf'], saved['weights'],
                                     float(saved['noise_std']), np.deg2rad(1), .5/g['field_A'], order=32)
    op.establish_coefficient_scaling()
    ordinary_forward, ordinary_adjoint = op._transform_forward, op._transform_adjoint
    start = time.perf_counter()
    forward = finufft.Plan(3, 3, n_trans=10, eps=op.eps, isign=1, nthreads=1)
    forward.setpts(*op.source, s=op.destination[0], t=op.destination[1], u=op.destination[2])
    adjoint = finufft.Plan(3, 3, n_trans=10, eps=op.eps, isign=1, nthreads=1)
    adjoint.setpts(*op.destination, s=op.source[0], t=op.source[1], u=op.source[2])
    setup = time.perf_counter()-start
    planned_forward = lambda strengths: forward.execute(np.ascontiguousarray(strengths))
    planned_adjoint = lambda strengths: adjoint.execute(np.ascontiguousarray(strengths))
    rng = np.random.default_rng(609871); records = []
    for iteration in range(6):
        v = rng.normal(size=op.shape[0]); v /= np.linalg.norm(v)
        # Include a weight update: plans depend on positions, not strengths.
        op.set_weights(saved['weights']*(1+.01*rng.normal(size=len(saved['weights']))))
        outputs = {}; times = {}
        for mode in (['ordinary', 'planned'] if iteration % 2 == 0 else ['planned', 'ordinary']):
            op._transform_forward = planned_forward if mode == 'planned' else ordinary_forward
            op._transform_adjoint = planned_adjoint if mode == 'planned' else ordinary_adjoint
            start = time.perf_counter(); outputs[mode] = op.spatial_gram(v); times[mode] = time.perf_counter()-start
        error = float(np.linalg.norm(outputs['planned']-outputs['ordinary'])/np.linalg.norm(outputs['ordinary']))
        records.append({'iteration': iteration, 'seconds': times, 'relative_discrepancy': error})
        print(records[-1], flush=True)
    result = {'complete': True, 'source_snapshot': source, 'particles': 128, 'radius': 5, 'order': 32,
              'plan_setup_seconds': setup, 'records': records,
              'median_ordinary_seconds': float(np.median([r['seconds']['ordinary'] for r in records])),
              'median_planned_seconds': float(np.median([r['seconds']['planned'] for r in records])),
              'maximum_relative_discrepancy': max(r['relative_discrepancy'] for r in records)}
    result['median_speedup'] = result['median_ordinary_seconds']/result['median_planned_seconds']
    path.write_text(json.dumps(result, indent=2)+'\n'); print('SPEEDUP', result['median_speedup'], flush=True)


if __name__ == '__main__': main()
