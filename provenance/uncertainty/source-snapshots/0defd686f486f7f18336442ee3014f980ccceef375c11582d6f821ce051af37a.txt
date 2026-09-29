#!/usr/bin/env python3
"""Continuous nonlinear-pose audit of saved fixed-pose weights (development)."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry, VoxelReference
from fourier_splats.uq_continuous import cell_forward, cell_target_coefficients, continuous_residual_norm
from fourier_splats.uq_continuous_pose import continuous_pose_audit, pose_cell_forward, perturbed_geometry
from fourier_splats.uq_provenance import source_snapshot
from fourier_splats.uncertainty import bias_aware_half_width
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
MAPS = {'10028': '2660', '10049': '6487', '10076': '8434'}
SNAPSHOT = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py'])


def run(args, dataset):
    source_path = BASE/args.source/f'{dataset}.json'
    source = json.loads(source_path.read_text()); config = source['config']
    g = particle_geometry(ROOT, dataset, 'inference_half0', radius=config['frequency_radius'],
                          count=config['particles'], seed=config['seed'])
    checkpoint = np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
    op, pilot, _, noise = model(g, checkpoint, 24)
    pilot_full = op.expand(pilot)
    nominal_pilot = cell_forward(g['k'], g['ctf'], pilot_full, 24, noise)
    reference = VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map', box=64).volume.ravel()
    reference /= np.linalg.norm(reference)
    out = BASE/args.output; out.mkdir(parents=True, exist_ok=True)
    result = {'stage': 'development: continuous L2 cube, bounded nonlinear pose, known fixed CTF and Gaussian whitening',
              'dataset': dataset, 'config': vars(args), 'source_snapshot': SNAPSHOT,
              'source_json_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'class': 'radius 2 around unit L2 constant-cell pilot; density need not be constant on cells',
              'weight_selection': 'saved continuous fixed-pose weights; nonlinear certificate is a post-audit, not an optimized pose solution',
              'records': []}
    rng = np.random.default_rng(args.seed)
    for name in args.targets.split(','):
        centers = [[0, 0, 0]] if name == 'center' else [[0, 0, .08], [0, 0, -.08]]
        signs = [1] if name == 'center' else [1, -1]
        saved_path = BASE/args.source/f'{dataset}-{name}-{args.width}-weights.npz'
        saved = np.load(saved_path); w = saved['weights']; np.testing.assert_array_equal(saved['indices'], g['indices'])
        np.testing.assert_allclose(saved['noise_std'], noise)
        fixed = next(r for r in source['targets'] if r['target'] == name and r['width_fraction_field'] == args.width)
        B = 2.; P = float(np.linalg.norm(pilot_full)); residual = fixed['fit']['bias']/B
        sd = float(np.linalg.norm(w)); no_data = B*fixed['fit']['target_norm']
        target_pilot = float(cell_target_coefficients(24, centers, signs, args.width)@pilot_full)
        true_reference = float(cell_target_coefficients(64, centers, signs, args.width)@reference)
        for angle_degrees in map(float, args.angles.split(',')):
            start = time.perf_counter(); angle = np.deg2rad(angle_degrees); shift = .01/24
            callback = lambda row: print(dataset, name, angle_degrees, row, flush=True)
            audit = continuous_pose_audit(g['k'], g['q'], g['ctf'], w, noise, angle, shift, B, P, residual,
                                           order=args.order, backend=args.backend, callback=callback)
            half = float(bias_aware_half_width(sd, audit['total_bias'])); choose_data = bool(half < no_data)
            row = {'target': name, 'width_fraction_field': args.width, 'angle_degrees': angle_degrees,
                   'shift_fraction_field': shift, 'noise_sd': sd, 'audit': audit, 'half_width_before_fallback': half,
                   'selected_half_width': min(half, no_data), 'selected_data_estimator': choose_data,
                   'width_fraction_of_no_data': min(half, no_data)/no_data,
                   'fixed_pose_width_fraction_of_no_data': fixed['width_fraction_of_no_data'],
                   'weights_sha256': hashlib.sha256(w.tobytes()).hexdigest(), 'checks': []}
            # Feasible stress draws are geometry-only and independent of image noise.
            poses = [np.zeros((len(g['k']), 5))]
            for _ in range(args.pose_draws):
                u = rng.normal(size=(len(g['k']), 5)); u /= np.linalg.norm(u, axis=1)[:, None]; poses.append(u)
            coherent = np.zeros_like(poses[0]); coherent[:, 0] = 1; poses.append(coherent)
            for index, pose in enumerate(poses):
                signal = pose_cell_forward(g['k'], g['q'], g['ctf'], reference, 64, noise, pose, angle, shift)
                actual_bias = float(target_pilot+w@(signal-nominal_pilot)-true_reference)
                coverage = float(norm.cdf((half-actual_bias)/sd)-norm.cdf((-half-actual_bias)/sd))
                check = {'pose_index': index, 'pose_kind': 'nominal' if index == 0 else ('coherent_x_boundary' if index == len(poses)-1 else 'independent_ball_boundary'),
                         'reference_bias': actual_bias, 'reference_coverage_before_fallback': coverage}
                if abs(actual_bias) > audit['total_bias']+1e-6:
                    raise AssertionError('Reference nonlinear bias exceeds continuous bound')
                # Expensive analytic sinc norm independently eliminates the FULL
                # continuous density ball at one nonzero, coherent allowed pose.
                if index == len(poses)-1 and not args.skip_exact_stress:
                    rotated, translations = perturbed_geometry(g['k'], g['q'], pose, angle, shift)
                    perturbed_ctf = g['ctf']*np.exp(2j*np.pi*translations)
                    exact = continuous_residual_norm(rotated, perturbed_ctf, w, noise, centers, signs, args.width)
                    pilot_signal = pose_cell_forward(g['k'], g['q'], g['ctf'], pilot_full, 24, noise, pose, angle, shift)
                    pilot_bias = float(w@(pilot_signal-nominal_pilot))
                    feasible_worst = abs(pilot_bias)+B*exact['residual_norm']
                    check.update(continuous_density_ball_eliminated_bias=feasible_worst,
                                 feasible_bias_fraction_of_upper=feasible_worst/audit['total_bias'])
                    if feasible_worst > audit['total_bias']+1e-6:
                        raise AssertionError('Continuous nonlinear feasible bias exceeds certificate')
                row['checks'].append(check)
            row['seconds'] = time.perf_counter()-start; result['records'].append(row)
            (out/f'{dataset}.json').write_text(json.dumps(result, indent=2)+'\n')
            print(dataset, name, angle_degrees, 'continuous pose width/no-data', row['width_fraction_of_no_data'],
                  'uncapped', half/no_data, 'seconds', row['seconds'], flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--datasets', default='10028,10049,10076')
    p.add_argument('--source', default='continuous-quadrature-optimized'); p.add_argument('--targets', default='center,contrast')
    p.add_argument('--width', type=float, default=.07); p.add_argument('--angles', default='.1,.5')
    p.add_argument('--order', type=int, default=32); p.add_argument('--backend', choices=['direct', 'nufft'], default='nufft')
    p.add_argument('--seed', type=int, default=609393); p.add_argument('--pose-draws', type=int, default=4)
    p.add_argument('--skip-exact-stress', action='store_true'); p.add_argument('--output', default='continuous-pose')
    args = p.parse_args()
    for dataset in args.datasets.split(','):
        run(args, dataset)
