#!/usr/bin/env python3
"""One actual-geometry operation-equivalence/timing check; no new interval."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_shift_only import ShiftAwarePolynomialPoseFieldOperator
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def main():
    path = BASE/'audit-regressions/shift-only-operation-equivalence.json'
    if path.exists(): raise RuntimeError('Preserve prior outcome')
    result = {'complete': False, 'scope': 'Operation-equivalence and single-call timing only; no spectral event or new inference',
              'source_snapshot': source_snapshot(ROOT, Path(__file__))}
    def save(): path.write_text(json.dumps(result, indent=2)+'\n')
    save()
    try:
        fit_path = BASE/'pilot-selected-fixed-v1/10049-pilot_region_1.json'; fit = json.loads(fit_path.read_text())
        row = fit['targets'][0]; wp = fit_path.with_name(f"10049-pilot_region_1-{row['width_fraction_field']}-weights.npz")
        saved = np.load(wp); g = particle_geometry(ROOT, '10049', 'inference_half0', radius=12, count=128, seed=fit['config']['seed'])
        np.testing.assert_array_equal(saved['indices'], g['indices'])
        kwargs = dict(order=80, nthreads=1)
        args = (g['k'], g['q'], g['ctf'], saved['weights'], float(saved['noise_std']), 0., .5/g['field_A'])
        original = PolynomialPoseFieldOperator(*args, **kwargs); fast = ShiftAwarePolynomialPoseFieldOperator(*args, **kwargs)
        scales = json.loads((BASE/'pilot-selected-pose-v1/10049-pilot_region_1-1.json').read_text())['pose_scaling']['column_denominators']
        original.denominators = np.array(scales); fast.denominators = np.array(scales)
        rng = np.random.default_rng(642001); u = rng.normal(size=original.shape[1]); v = rng.normal(size=original.shape[0])
        times = {}; values = {}
        for label, op in [('original', original), ('reduced', fast)]:
            start = time.perf_counter(); forward = op.matvec(u); times[label+'_forward_seconds'] = time.perf_counter()-start
            start = time.perf_counter(); adjoint = op.rmatvec(v); times[label+'_adjoint_seconds'] = time.perf_counter()-start
            values[label] = (forward, adjoint)
        checks = {}
        for index, name in enumerate(['forward', 'adjoint']):
            old, new = values['original'][index], values['reduced'][index]
            relative = float(np.linalg.norm(old-new)/max(np.linalg.norm(old), 1e-300))
            checks[name] = {'relative_norm_difference': relative, 'maximum_absolute_difference': float(np.max(abs(old-new))),
                            'passed': bool(relative < 1e-9)}
        result.update(times=times, comparisons=checks, dataset='10049', particles=128, frequencies_per_particle=original.nq,
            quadrature_order=80, shape=original.shape, threads=1, rotation_radius=0., shift_radius_A=.5,
            seed=642001, fit_sha256=hashlib.sha256(fit_path.read_bytes()).hexdigest(),
            weights_sha256=hashlib.sha256(wp.read_bytes()).hexdigest(),
            speedup=(times['original_forward_seconds']+times['original_adjoint_seconds'])/(times['reduced_forward_seconds']+times['reduced_adjoint_seconds']),
            timing_limit='Single call per direction, concurrent workloads, no steady-state or pipeline throughput claim')
        if not all(r['passed'] for r in checks.values()): raise AssertionError('Actual geometry operation mismatch')
        result['complete'] = True; save(); print(json.dumps({k:result[k] for k in ['times','comparisons','speedup']}), flush=True)
    except Exception as exc: result.update(error=repr(exc)); save(); raise


if __name__ == '__main__': main()
