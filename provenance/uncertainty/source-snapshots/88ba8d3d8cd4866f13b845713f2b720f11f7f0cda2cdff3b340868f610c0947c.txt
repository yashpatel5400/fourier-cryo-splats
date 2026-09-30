#!/usr/bin/env python3
"""C5: exact-continuous selected-block trace diagnostics; no scale retuning."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from fourier_splats.uq_cubic_pose import CubicPoseFieldOperator, SPATIAL
from fourier_splats.uq_polynomial_trace import polynomial_fourier_block_trace
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'
POSITIONS = np.array([0, 64, 127])


def main():
    path = BASE/'audit-regressions/cubic-design-scale-check.json'
    if path.exists(): raise RuntimeError('Preserve previous outcome')
    start = time.perf_counter()
    result = {'complete': False, 'scope': 'Three deterministic particles, analytic continuous block traces versus design scales; no scale or interval is changed',
              'positions': POSITIONS.tolist(), 'source_snapshot': source_snapshot(ROOT, Path(__file__))}
    def save(): path.write_text(json.dumps(result, indent=2)+'\n')
    save()
    try:
        audit_path = BASE/'cubic-pose-probe/10049-pilot_region_1-1.json'; audit = json.loads(audit_path.read_text())
        fit_path = BASE/'pilot-selected-fixed-v1/10049-pilot_region_1.json'; fit = json.loads(fit_path.read_text()); row = fit['targets'][0]
        wp = fit_path.with_name(f"10049-pilot_region_1-{row['width_fraction_field']}-weights.npz")
        if hashlib.sha256(fit_path.read_bytes()).hexdigest() != audit['source_fit_sha256'] or hashlib.sha256(wp.read_bytes()).hexdigest() != audit['source_weights_sha256']:
            raise AssertionError('Changed estimator')
        g = particle_geometry(ROOT, '10049', 'inference_half0', radius=12, count=128, seed=fit['config']['seed'])
        saved = np.load(wp); np.testing.assert_array_equal(g['indices'], saved['indices']); n, nq = g['ctf'].shape
        op = CubicPoseFieldOperator(g['k'][POSITIONS], g['q'][POSITIONS], g['ctf'][POSITIONS],
            saved['weights'].reshape(n, 2*nq)[POSITIONS], float(saved['noise_std']), np.deg2rad(1.), .5/g['field_A'], order=12, nthreads=1)
        reproduced = op.establish_quadrature_design_scaling(12)['group_scales'].reshape(3, 3)
        stored = np.array(audit['pose_scaling']['group_scales']).reshape(3, 128)[:, POSITIONS]
        np.testing.assert_allclose(reproduced, stored, rtol=1e-12)
        result.update(source_fit_sha256=audit['source_fit_sha256'], source_weights_sha256=audit['source_weights_sha256'], records=[])
        save()
        for i, position in enumerate(POSITIONS):
            coefficients = op.polynomials[i]*op.c[i, :, None, None]
            for degree, (a, b) in enumerate(((0, 5), (5, 20), (20, 55)), 1):
                record = polynomial_fourier_block_trace(op.k[i], coefficients[:, a:b], SPATIAL)
                scale = float(stored[degree-1, i])
                result['records'].append(dict(particle_position=int(position), degree=degree, design_scale=scale,
                    continuous_norm_over_design_scale=record['norm_unpadded']/scale, **record))
                save(); print(position, degree, result['records'][-1]['continuous_norm_over_design_scale'], flush=True)
        result.update(complete=True, seconds=time.perf_counter()-start,
            interpretation='Trace ratios measure selected scaling efficiency; neither prove optimal block scales nor bound all floating-point error')
        save()
    except Exception as exc: result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__': main()
