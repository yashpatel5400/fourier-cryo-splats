#!/usr/bin/env python3
"""Post-review-candidate development: select bandwidth using a pose audit.

No frozen experiment is modified. Geometry, target, pilot and physical noise
scale stay fixed across bandwidths. Selection uses only design-derived width;
reference values are diagnostic outcomes, never selection criteria.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry, VoxelReference
from fourier_splats.uq_continuous import (continuous_certificate, continuous_residual_norm,
                                        cell_forward, cell_target_coefficients)
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_continuous_pose import continuous_pose_audit, pose_cell_forward
from fourier_splats.uq_continuous_moments import integrated_cubic_remainder
from fourier_splats.uncertainty import bias_aware_half_width
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
MAPS = {'10028': '2660', '10049': '6487', '10076': '8434'}


def run(args, dataset, snapshot):
    path = BASE/args.output/f'{dataset}.json'
    if path.exists():
        raise RuntimeError('Preserve prior outcomes; choose another output name')
    path.parent.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    full = particle_geometry(ROOT, dataset, 'inference_half0', radius=5, count=128, seed=609315)
    checkpoint = np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
    op, pilot, _, noise = model(full, checkpoint, 24)
    pilot = op.expand(pilot); P = float(np.linalg.norm(pilot)); B = 2.
    rho = VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map', box=64).volume.ravel()
    rho /= np.linalg.norm(rho)
    result = {'stage': 'post-candidate exploratory bandwidth choice; not a frozen confirmation',
              'dataset': dataset, 'config': vars(args), 'source_snapshot': snapshot,
              'physical_noise_std_fixed_across_bands': float(noise), 'density_radius': B,
              'pilot_norm': P, 'geometry_indices': full['indices'].tolist(),
              'geometry_sha256': hashlib.sha256(full['k'].tobytes()+full['ctf'].tobytes()).hexdigest(),
              'reference_sha256': hashlib.sha256(rho.tobytes()).hexdigest(),
              'selection_rule': 'Minimum audited width over the declared band grid and no-data estimator; no reference or inference noise used.',
              'records': [], 'complete': False}
    path.write_text(json.dumps(result, indent=2)+'\n')
    for radius in map(float, args.radii.split(',')):
        keep = np.linalg.norm(full['q'][0], axis=1) <= radius+1e-10
        if radius > 5 or not keep.any():
            raise ValueError('Probe bands must be nonempty subsets of radius five')
        g = {key: full[key][:, keep] for key in ['k', 'q', 'ctf']}
        gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=40, preconditioner_rank=1024)
        nominal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
        for target in args.targets.split(','):
            centers = [[0, 0, 0]] if target == 'center' else [[0, 0, .08], [0, 0, -.08]]
            signs = [1] if target == 'center' else [1, -1]
            fit = continuous_certificate(gram, centers, signs, .07, B, maxiter=150, rtol=.005)
            w = fit.pop('weights'); sd = float(np.linalg.norm(w))
            exact = continuous_residual_norm(g['k'], g['ctf'], w, noise, centers, signs, .07)
            if abs(exact['residual_norm']-fit['bias']/B) > 1e-5:
                raise AssertionError('Independent continuous integral check failed')
            pilot_target = float(cell_target_coefficients(24, centers, signs, .07)@pilot)
            truth = float(cell_target_coefficients(64, centers, signs, .07)@rho)
            np.savez(path.parent/f'{dataset}-{target}-radius{radius:g}-weights.npz', weights=w,
                     indices=full['indices'], noise_std=noise, frequency_mask=keep)
            for degrees in map(float, args.angles.split(',')):
                tick = time.perf_counter(); angle = np.deg2rad(degrees); shift = .01/24
                audit = continuous_pose_audit(g['k'], g['q'], g['ctf'], w, noise, angle, shift,
                                              B, P, exact['residual_norm'], order=32)
                moment = integrated_cubic_remainder(g['k'], g['q'], g['ctf'], w, noise, angle, shift, B, P)
                moment['per_particle_operator_third_bound'] = moment['per_particle_operator_third_bound'].tolist()
                bias_bound = audit['density_bias']+audit['pose_polynomial_bias']+min(audit['remainder_bias'], moment['remainder_bias'])
                half = float(bias_aware_half_width(sd, bias_bound)); no_data = B*fit['target_norm']
                use_no_data = half >= no_data; selected_half = min(half, no_data)
                checks = []
                for scenario in ['nominal', 'coherent_x', 'random_boundary']:
                    pose = np.zeros((len(g['k']), 5))
                    if scenario == 'coherent_x':
                        pose[:, 0] = 1.
                    elif scenario == 'random_boundary':
                        pose = np.random.default_rng(609571+int(dataset)).normal(size=pose.shape)
                        pose /= np.linalg.norm(pose, axis=1)[:, None]
                    signal = pose_cell_forward(g['k'], g['q'], g['ctf'], rho, 64, noise, pose, angle, shift)
                    mean = pilot_target if use_no_data else float(pilot_target+w@(signal-nominal))
                    actual_sd = 0. if use_no_data else sd; actual_bias = mean-truth
                    coverage = (float(norm.cdf((selected_half-actual_bias)/actual_sd)-norm.cdf((-selected_half-actual_bias)/actual_sd))
                                if actual_sd else float(abs(actual_bias) <= selected_half))
                    power = (float(norm.cdf((np.sign(truth)*mean-selected_half)/actual_sd))
                             if actual_sd else float(np.sign(truth)*mean > selected_half))
                    checks.append({'scenario': scenario, 'true_target': truth, 'expected_center': mean,
                                   'actual_bias': actual_bias, 'analytic_coverage': coverage,
                                   'correct_sign_probability': power})
                row = {'target': target, 'frequency_radius': radius, 'angle_degrees': degrees,
                       'fit': fit, 'audit': audit, 'moment_remainder': moment,
                       'total_bias_bound': bias_bound, 'half_width_before_fallback': half,
                       'relative_half_width': selected_half/no_data, 'uses_no_data': use_no_data,
                       'reference_checks': checks, 'audit_seconds': time.perf_counter()-tick}
                result['records'].append(row)
                path.write_text(json.dumps(result, indent=2)+'\n')
                print(dataset, target, radius, degrees, 'relative width', row['relative_half_width'],
                      'min reference power', min(c['correct_sign_probability'] for c in checks), flush=True)
        del gram
    result['selected'] = []
    for target in args.targets.split(','):
        for degrees in map(float, args.angles.split(',')):
            rows = [r for r in result['records'] if r['target']==target and r['angle_degrees']==degrees]
            best = min(rows, key=lambda r: r['relative_half_width'])
            result['selected'].append({'target': target, 'angle_degrees': degrees,
                                      'frequency_radius': best['frequency_radius'],
                                      'relative_half_width': best['relative_half_width'],
                                      'minimum_reference_power': min(c['correct_sign_probability'] for c in best['reference_checks'])})
    result['complete'] = True; result['seconds'] = time.perf_counter()-start
    path.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--datasets', default='10028')
    parser.add_argument('--targets', default='center')
    parser.add_argument('--radii', default='2,3,4,5')
    parser.add_argument('--angles', default='.5,1,2')
    parser.add_argument('--output', default='pose-aware-bandwidth-probe')
    args = parser.parse_args()
    snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py'])
    for dataset in args.datasets.split(','):
        run(args, dataset, snapshot)
