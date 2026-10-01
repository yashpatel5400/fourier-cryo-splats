#!/usr/bin/env python3
"""Numerical residual energy for a saved nonconverged Gaussian fit."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_provenance import source_snapshot


def main():
    root = Path(__file__).resolve().parents[1]; base = root/'results/uncertainty/development'
    out = base/'audit-regressions/continuous-gaussian-solve-gap-first-case.json'
    if out.exists():
        raise RuntimeError('Preserve previous diagnostic')
    source = base/'continuous-gaussian-review2-v1/10028.json'
    data = json.loads(source.read_text()); row = data['records'][0]
    if row['target'] != 'pilot_region_1' or row['pose_model'] != 'fixed_pose' or row['prior_directional_sd'] != 2.:
        raise ValueError('Unexpected first fit')
    path = base/'continuous-gaussian-review2-v1/10028-pilot_region_1-tau2-fixed_pose.npz'
    if hashlib.sha256(path.read_bytes()).hexdigest() != row['weights_sha256']:
        raise ValueError('Saved fit changed')
    saved = np.load(path)
    g = particle_geometry(root, '10028', 'inference_half0', radius=12, count=128, seed=609315)
    np.testing.assert_array_equal(g['indices'], saved['indices'])
    tick = time.perf_counter()
    gram = QuadratureObservationGram(g['k'], g['ctf'], float(saved['noise_std']), order=80, preconditioner_rank=1024)
    gram.nthreads = 1
    a, ell2 = gram.target(row['centers'], [1], row['width_fraction_field'])
    w = saved['weights']; tau = 2.; residual = gram.matvec(w)+w/tau**2-a
    energy = float(residual@(gram.preconditioner(tau**-2)@residual))
    if energy < 0:
        raise ArithmeticError('Negative residual energy')
    integration_error = gram.quadrature_error(w)['gram_action_norm']
    gap = tau**2*(np.sqrt(energy)+tau*integration_error)**2
    result = dict(complete=True, seconds=time.perf_counter()-tick, saved_fit_record=row,
        weights_sha256=row['weights_sha256'],
        source_snapshot=source_snapshot(root, Path(__file__), ['research/uncertainty/GAUSSIAN-SOLVE-GAP-NOTE.md']),
        relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(a)),
        residual_preconditioned_energy=energy, quadrature_residual_norm_error=integration_error,
        variance_gap_upper_diagnostic=float(gap),
        trivial_variance_gap_upper_diagnostic=float(tau**4*(np.linalg.norm(residual)+integration_error)**2),
        relative_variance_gap_upper_diagnostic=float(gap/row['fit']['variational_variance_upper']),
        scope='Real-arithmetic residual quadratic identity with analytic quadrature control. Floating-point factor/NUFFT errors unenclosed. Original CG failure unchanged; this does not establish an exact posterior mean.')
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['relative_residual'], result['relative_variance_gap_upper_diagnostic'], flush=True)


if __name__ == '__main__':
    main()
