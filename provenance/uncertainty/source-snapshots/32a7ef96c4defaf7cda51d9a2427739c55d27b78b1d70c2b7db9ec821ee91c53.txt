#!/usr/bin/env python3
"""One declared fixed-weight cubic-pose audit; no adaptive empirical cases."""
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
import numpy as np
import finufft
from fourier_splats.uq_cubic_pose import (CubicPoseFieldOperator, cell_moments_cubic,
    direct_cell_moments_cubic, cubic_residual_cross, BLOCK_SIZES)
from fourier_splats.uq_higher_remainder import higher_order_particle_remainders
from fourier_splats.uq_continuous_pose import polynomial_kernel_error
from fourier_splats.uq_joint_bias import joint_density_pose_bias
from fourier_splats.uq_continuous import cell_target_coefficients
from fourier_splats.uq_random_spectral import gaussian_power_upper
from fourier_splats.uq_intervals import bias_aware_half_width_stable, reference_interval_summary
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model
from probe_uq_ball_remainder import validate_source_class

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
SOURCE = BASE/'pilot-selected-pose-v1/10049-pilot_region_1-1.json'
CONFIG = {'order': 80, 'design_order': 12, 'iterations': 40, 'probes': 4,
          'certificate_seed': 640001, 'delta': 1e-6/12, 'threads': 2, 'degree': 3}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if 'tmp/nufft-openmp' not in finufft.__file__ or 'torch' in sys.modules:
        raise RuntimeError('Isolated CPU-only threaded FINUFFT required')
    out = BASE/'cubic-pose-probe'; out.mkdir(parents=True, exist_ok=True)
    path = out/SOURCE.name
    if path.exists():
        raise RuntimeError('Preserve previous outcome')
    old = json.loads(SOURCE.read_text()); fit_path = ROOT/old['config']['fit']
    fit = json.loads(fit_path.read_text())
    if not old.get('complete') or not fit.get('complete') or old.get('error') or fit.get('error'):
        raise ValueError('Completed source audit and fit required')
    B, P, origin = validate_source_class(old, fit)
    if (B, P) != (2., 1.) or old['dataset'] != '10049' or old['target'] != 'pilot_region_1':
        raise ValueError('Changed declared case')
    if old['config']['angle'] != 1. or old['config']['shift_A'] != .5:
        raise ValueError('Changed pose radius')
    cfg = fit['config']; width = old['width_fraction_field']; row = fit['targets'][0]
    if cfg['particles'] != 128 or cfg['frequency_radius'] != 12 or len(fit['targets']) != 1:
        raise ValueError('Changed declared cohort, frequency band or target')
    if sha(fit_path) != old['source_fit_sha256']:
        raise ValueError('Source fit changed')
    wp = fit_path.with_name(f'10049-pilot_region_1-{width}-weights.npz')
    if sha(wp) != old['source_weights_sha256']:
        raise ValueError('Source weights changed')
    result = {'complete': False, 'config': CONFIG, 'dataset': '10049', 'target': 'pilot_region_1',
        'scope': 'Prespecified single-case development diagnostic, conditional known-noise cubic pose audit; not experimental calibration or pose-optimized weights',
        'source_audit': str(SOURCE.relative_to(ROOT)), 'source_audit_sha256': sha(SOURCE),
        'source_fit_sha256': sha(fit_path), 'source_weights_sha256': sha(wp),
        'density_radius': B, 'pilot_norm_bound': P, 'class_metadata_origin': origin,
        'pose_set': old['pose_set'], 'alpha_total': old['alpha_total'], 'alpha_numerical': CONFIG['delta'],
        'alpha_noise': old['alpha_total']-CONFIG['delta'],
        'numerical_scope': 'Ordinary floating-point guards and independent numerical checks, not interval arithmetic',
        'source_snapshot': source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py',
            'scripts/probe_uq_ball_remainder.py', 'research/uncertainty/CUBIC-POSE-PROBE.md',
            'research/uncertainty/HIGHER-ORDER-REMAINDER-PROBE.md'])}
    start = time.perf_counter()
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    def progress(record):
        print(record, flush=True)
    save()
    try:
        g = particle_geometry(ROOT, '10049', 'inference_half0', radius=12, count=128, seed=cfg['seed'])
        saved = np.load(wp); np.testing.assert_array_equal(saved['indices'], g['indices'])
        w = saved['weights']; noise = float(saved['noise_std']); angle = np.deg2rad(1.); shift = .5/g['field_A']
        op = CubicPoseFieldOperator(g['k'], g['q'], g['ctf'], w, noise, angle, shift,
            order=80, nthreads=2, block_particles=16)
        scaling = op.establish_quadrature_design_scaling(12, callback=progress)
        scales = np.asarray(scaling['group_scales'], dtype='<f8'); den = np.asarray(op.denominators, dtype='<f8')
        assert scales.shape == (3*128,) and den.shape == (55*128,) and (scales > 0).all()
        result.update(particles=128, fourier_pairs_per_particle=op.nq, target_width_A=width*g['field_A'],
            coefficient_tensor_bytes=op.polynomials.nbytes,
            pose_scaling={'block_sizes': BLOCK_SIZES, 'group_scales': scales.tolist(),
                'column_denominators': den.tolist(), 'group_scales_sha256': hashlib.sha256(scales.tobytes()).hexdigest(),
                'column_denominators_sha256': hashlib.sha256(den.tobytes()).hexdigest(),
                'encoding': 'C-order little-endian IEEE-754 float64'}, scaling_method=scaling['scaling_method'])
        # Recover raw amplitudes directly, including zero-CTF frequencies.
        wr = w.reshape(op.n, 2*op.nq); amp = np.hypot(wr[:, :op.nq], wr[:, op.nq:])
        masses = np.einsum('naj,nj->na', op.coefficient_bound_matrix(), amp)
        integration_pad = float(polynomial_kernel_error(g['k'], 80, degree=6)*np.sum(masses*masses))
        result.update(setup_seconds=time.perf_counter()-start, integration_eigenvalue_pad=integration_pad)
        save()
        spectral = gaussian_power_upper(op.spatial_gram, op.shape[0], CONFIG['delta'], 4, 40, 640001, callback=progress)
        L = float(np.sqrt(scales.sum())); field = L*np.sqrt(spectral['eigenvalue_upper']+integration_pad)
        result['spectral_upper_certificate'] = spectral
        remainders = [higher_order_particle_remainders(g['k'][i], g['q'][i], op.c[i], angle, shift,
                       degrees=(3,), domain='cube') for i in range(op.n)]
        remainder_field = np.array([r['records'][0]['field_remainder'] for r in remainders])
        cross = cubic_residual_cross(op, w, noise, old['centers_fraction_field'], old['signs'], width)
        cross_vector = cross.pop('cross_vector')
        checkpoint = np.load(BASE/'representation/10049/real_particles-spacing-2.0.npz')
        cellop, pilot, _, check_noise = model(g, checkpoint, 24, noise=noise); pilot = cellop.expand(pilot)
        np.testing.assert_allclose(check_noise, noise); np.testing.assert_allclose(np.linalg.norm(pilot), P, rtol=1e-12)
        moments = cell_moments_cubic(g['k'], pilot, 24, nthreads=2)
        pilot_vector = op.pair_moments(moments)
        positions = np.array([0, 64, 127]); direct = direct_cell_moments_cubic(g['k'][positions], pilot, 24)
        relative = float(np.linalg.norm(direct-moments[positions])/max(np.linalg.norm(direct), 1e-300))
        absolute = float(np.max(abs(direct-moments[positions])))
        result['pilot_numerical_check'] = {'particle_positions': positions.tolist(),
            'relative_moment_difference': relative, 'maximum_absolute_moment_difference': absolute,
            'passed': bool(relative < 1e-8 or absolute < 1e-12),
            'scope': 'Selected frequencies: independent physical-cell sums and 24-node within-cell quadrature; not global numerical certification'}
        if not result['pilot_numerical_check']['passed']:
            raise AssertionError('Independent pilot moment check failed')
        joint = joint_density_pose_bias(row['fit']['bias']/B, field, cross['norm_upper'], L, B, P,
                                        float(remainder_field.sum()), float(np.linalg.norm(pilot_vector)))
        half = bias_aware_half_width_stable(old['noise_sd'], joint['bias_upper'], result['alpha_noise'])
        no_data = B*row['fit']['target_norm']
        pilot_target = float(cell_target_coefficients(24, old['centers_fraction_field'], old['signs'], width)@pilot)
        np.testing.assert_allclose(pilot_target, old['pilot_target'], rtol=1e-12)
        checks = [dict(scenario=r['scenario'], raw_expected_center=r['raw_expected_center'],
                       **reference_interval_summary(r['true_target'], r['raw_expected_center'], old['noise_sd'], pilot_target, half, no_data))
                  for r in old['reference_checks']]
        np.savez(path.with_suffix('.npz'), residual_pose_cross=cross_vector, pilot_pose_cross=pilot_vector,
            quartic_remainder_per_particle=remainder_field, sobolev_norms=np.array([r['sobolev_norms'] for r in remainders]),
            group_scales=scales, column_denominators=den, indices=g['indices'])
        result.update(complete=True, bound=joint, cross=cross, lifted_radius=L, pose_field_upper=float(field),
            density_bias=row['fit']['bias'], remainder_bias=(B+P)*float(remainder_field.sum()),
            noise_sd=old['noise_sd'], half_width=half, selected_relative_half_width=min(1., half/no_data),
            uses_no_data=bool(half >= no_data), reference_checks=checks, pilot_target=pilot_target,
            reference_checks_origin='Exact nonlinear source means reused at identical weights, geometry, noise and target; no new reference projections or outcomes selected',
            minimum_reference_power=min(r['correct_sign_probability'] for r in checks),
            reference_coverage_failure=any(r['analytic_coverage'] < 1-result['alpha_noise']-1e-8 for r in checks),
            arrays_sha256=sha(path.with_suffix('.npz')), seconds=time.perf_counter()-start,
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024))
        save(); print('COMPLETE', result['selected_relative_half_width'], result['minimum_reference_power'], flush=True)
    except Exception as exc:
        result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
