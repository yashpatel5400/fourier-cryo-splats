#!/usr/bin/env python3
"""Same-weight exploratory enclosing-ball remainder probe; retain all outcomes."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from fourier_splats.uq_ball_remainder import particle_ball_remainder
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_joint_bias import sharp_cube_cubic_coefficients, joint_density_pose_bias, scaled_pose_radius
from fourier_splats.uq_intervals import bias_aware_half_width_stable, reference_interval_summary
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_source_class(audit, fit_record):
    """Fail closed on a different/unknown class; preserve explicit legacy metadata."""
    expected = 'Per-particle joint ball: squared scaled rotation norm plus squared translation norm <= 1'
    legacy_pose = 'Per-particle joint ball: squared scaled rotation norm plus squared scaled translation norm <= 1'
    if audit.get('pose_set') not in {expected, legacy_pose}:
        raise ValueError('This probe requires an explicitly recorded joint pose ball')
    if 'density_radius' in audit and 'pilot_norm_bound' in audit:
        B, P = float(audit['density_radius']), float(audit['pilot_norm_bound'])
        origin = 'Explicit numeric source-audit fields'
    elif fit_record.get('class') == 'radius 2 around a unit-L2 constant-cell pilot; independent unit-L2 64-cell reference generator':
        B, P = 2., 1.
        origin = 'Exact legacy source-fit class string; source audit has no numeric radius fields'
    else:
        raise ValueError('Source density class is missing or unknown')
    if not np.isfinite([B, P]).all() or B <= 0 or P < 0:
        raise ValueError('Invalid source density radii')
    return B, P, origin


def recompute_bias(audit, B, P, cubic):
    """Re-evaluate both branches, rather than assume an additive decomposition."""
    L = scaled_pose_radius(np.asarray(audit['pose_scaling']['group_scales']))
    return joint_density_pose_bias(audit['density_bias']/B, audit['pose_polynomial_bias']/(B+P),
        audit['cross']['norm_upper'], L, B, P, cubic/(B+P), audit['bound']['pilot_cross_upper'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--audit', default='results/uncertainty/development/continuous-high-band-pose-1024-10A/10049-center-1.json')
    parser.add_argument('--output', default='ball-remainder-probe')
    parser.add_argument('--domain', choices=['ball', 'cube'], default='ball')
    args = parser.parse_args()
    audit_path = ROOT/args.audit; audit = json.loads(audit_path.read_text())
    if not audit.get('complete') or audit.get('error') or not audit['config'].get('pilot_pairing'):
        raise ValueError('Completed pilot-refined high-band audit required')
    fit_path = ROOT/audit['config']['fit']; fit_record = json.loads(fit_path.read_text())
    if sha(fit_path) != audit['source_fit_sha256']:
        raise AssertionError('Source fit changed')
    B, P, class_origin = validate_source_class(audit, fit_record)
    reconstructed = recompute_bias(audit, B, P, audit['cubic_bias'])
    for key in ['bias_upper', 'joint_bias_upper', 'triangle_bias_upper', 'pilot_polynomial_bias_upper']:
        np.testing.assert_allclose(reconstructed[key], audit['bound'][key], rtol=2e-12, atol=1e-12)
    dataset = audit['dataset']; target = audit['target']; width = audit['width_fraction_field']
    weights_path = fit_path.with_name(f'{dataset}-{target}-{width}-weights.npz')
    if sha(weights_path) != audit['source_weights_sha256']:
        raise AssertionError('Source weights changed')
    out = ROOT/'results/uncertainty/development'/args.output; out.mkdir(parents=True, exist_ok=True)
    path = out/audit_path.name
    if path.exists():
        raise RuntimeError('Preserve previous outcome')
    start = time.perf_counter()
    result = {'complete': False, 'scope': 'Exploratory same-weight conditional pose remainder refinement, not experimental calibration',
              'config': vars(args), 'dataset': dataset, 'target': target,
              'source_audit_sha256': sha(audit_path), 'source_fit_sha256': sha(fit_path),
              'source_weights_sha256': sha(weights_path),
              'density_radius': B, 'pilot_norm_bound': P, 'class_metadata_origin': class_origin,
              'joint_ball': True, 'pose_set': audit['pose_set'], 'domain': args.domain,
              'recomputed_source_bound': reconstructed,
              'numerical_scope': 'Analytic identities plus heuristic magnitude guards; not interval arithmetic',
              'source_snapshot': source_snapshot(ROOT, Path(__file__), ['research/uncertainty/BALL-SOBOLEV-REMAINDER-PROPOSAL.md'])}
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    save()
    try:
        cfg = fit_record['config']
        g = particle_geometry(ROOT, dataset, 'inference_half0', radius=cfg['frequency_radius'], count=cfg['particles'], seed=cfg['seed'])
        saved = np.load(weights_path); np.testing.assert_array_equal(g['indices'], saved['indices'])
        noise = float(saved['noise_std']); n, nq = g['ctf'].shape
        w = saved['weights'].reshape(n, 2*nq)
        c = (w[:, :nq]+1j*w[:, nq:])*g['ctf']/noise
        angle = np.deg2rad(audit['config']['angle']); shift = audit['config']['shift_A']/g['field_A']
        old = np.sum(sharp_cube_cubic_coefficients(g['k'], g['q'], g['ctf']/noise, angle, shift, B+P)*np.hypot(w[:, :nq], w[:, nq:]), axis=1)
        np.testing.assert_allclose(sum(old), audit['cubic_bias'], rtol=1e-12)
        records = []
        for i in range(n):
            records.append(particle_ball_remainder(g['k'][i], g['q'][i], c[i], angle, shift, joint_ball=True, domain=args.domain))
            if (i+1) % 128 == 0:
                print('particles', i+1, 'seconds', time.perf_counter()-start, flush=True)
        ball = (B+P)*np.array([r['field_remainder'] for r in records])
        selected = np.minimum(old, ball)
        new_bound = recompute_bias(audit, B, P, float(sum(selected)))
        bias = new_bound['bias_upper']
        half = bias_aware_half_width_stable(audit['noise_sd'], bias, audit['alpha_noise'])
        row = next(r for r in fit_record['targets'] if r['target'] == target)
        no_data = B*row['fit']['target_norm']
        # All source scenarios used the no-data fallback; their center records the pilot target.
        if not audit['uses_no_data']:
            raise ValueError('This probe requires a source fallback record to recover its pilot target')
        checks = [dict(scenario=r['scenario'], **reference_interval_summary(r['true_target'], r['raw_expected_center'], audit['noise_sd'], r['expected_center'], half, no_data)) for r in audit['reference_checks']]
        np.savez(path.with_suffix('.npz'), old_cubic_per_particle=old, ball_cubic_per_particle=ball,
                 selected_cubic_per_particle=selected, sobolev_norms=np.array([r['sobolev_norms'] for r in records]),
                 path_speed_bounds=np.array([r['path_speed_bound'] for r in records]), indices=g['indices'])
        result.update(complete=True, seconds=time.perf_counter()-start, particles=n,
                      domain=args.domain,
                      old_cubic_bias=float(sum(old)), ball_cubic_bias=float(sum(ball)), selected_cubic_bias=float(sum(selected)),
                      enclosing_cubic_bias=float(sum(ball)), bound=new_bound,
                      particles_improved=int(np.sum(ball < old)), minimum_ball_over_old=float(np.min(ball/old)),
                      maximum_ball_over_old=float(np.max(ball/old)), bias_upper=bias, half_width=half,
                      selected_relative_half_width=min(1., half/no_data), uses_no_data=bool(half >= no_data),
                      reference_checks=checks, minimum_reference_power=min(r['correct_sign_probability'] for r in checks),
                      maximum_embedding_phase_residual=max(r['maximum_phase_residual'] for r in records),
                      embedding_residual_bias=(B+P)*sum(r['embedding_residual_remainder'] for r in records),
                      maximum_relative_sobolev_pad=max(d['roundoff_pad']/max(abs(d['squared_unpadded']), np.finfo(float).tiny) for r in records for d in r['sobolev_diagnostics']),
                      arrays_sha256=sha(path.with_suffix('.npz')))
        save(); print('DONE', result['old_cubic_bias'], result['selected_cubic_bias'], result['selected_relative_half_width'], flush=True)
    except Exception as exc:
        result.update(complete=False, error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
