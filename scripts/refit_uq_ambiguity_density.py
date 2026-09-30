#!/usr/bin/env python3
"""Return to the continuous density class after a feasible local pose step."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous import cell_forward
from fourier_splats.uq_continuous_pose import perturbed_geometry
from fourier_splats.uq_planned_transforms import PlannedQuadratureObservationGram
from fourier_splats.uq_two_pose_modulus import PhaseRotatedGram, optimize_two_pose_modulus
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(); p.add_argument('--pose-source', required=True)
    p.add_argument('--output', required=True); p.add_argument('--maxiter', type=int, default=25); args = p.parse_args()
    pose_path = BASE/args.pose_source; pose_record = json.loads(pose_path.read_text())
    if not pose_record.get('complete') or pose_record.get('error'):
        raise ValueError('Completed feasible pose source required')
    source = BASE/pose_record['config']['source']; prior = json.loads(source.read_text())
    if sha(source) != pose_record['source_record_sha256'] or sha(source.with_suffix('.npz')) != pose_record['source_arrays_sha256']:
        raise ValueError('Continuous source changed')
    if sha(pose_path.with_suffix('.npz')) != pose_record['arrays_sha256']:
        raise ValueError('Pose source changed')
    if prior['pose_set'] != 'joint_scaled_five_dimensional_ball':
        raise ValueError('Joint pose source required')
    pose_arrays = np.load(pose_path.with_suffix('.npz')); poses = pose_arrays['poses']
    maximum_pose = float(np.linalg.norm(poses, axis=-1).max())
    if maximum_pose > 1.:
        raise ValueError('Pose source outside class')
    out = BASE/args.output; out.mkdir(parents=True, exist_ok=True); path = out/source.name
    if path.exists():
        raise RuntimeError('Preserve prior refit')
    start = time.perf_counter()
    fields = ['dataset', 'target', 'rotation_radius_degrees', 'translation_radius_A', 'pose_set', 'density_radius',
              'pilot_norm', 'pilot_checkpoint_sha256', 'noise_std', 'field_A', 'width_fraction_field', 'alpha',
              'distance_budget', 'source_geometry_config']
    result = {k: prior[k] for k in fields}
    result.update(complete=False, config=vars(args), stage='Continuous density refit at locally optimized feasible poses; not global pose optimum',
                  source_pose_sha256=sha(pose_path), source_continuous_sha256=sha(source), maximum_scaled_pose_norm=maximum_pose,
                  source_snapshot=source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py', 'research/uncertainty/POSE-OPTIMIZED-AMBIGUITY-PROTOCOL.md']))
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    save()
    try:
        dataset = prior['dataset']; target = prior['target']; noise = prior['noise_std']
        g = particle_geometry(ROOT, dataset, 'inference_half0', radius=5, count=128, seed=prior['source_geometry_config']['seed'])
        initial = np.load(source.with_suffix('.npz'))
        np.testing.assert_array_equal(initial['indices'], g['indices']); np.testing.assert_array_equal(pose_arrays['indices'], g['indices'])
        checkpoint = BASE/'representation'/dataset/'real_particles-spacing-2.0.npz'
        if sha(checkpoint) != prior['pilot_checkpoint_sha256']:
            raise ValueError('Pilot changed')
        cellop, pilot, _, check_noise = model(g, np.load(checkpoint), 24); pilot = cellop.expand(pilot)
        np.testing.assert_allclose(noise, check_noise); np.testing.assert_allclose(np.linalg.norm(pilot), prior['pilot_norm'], rtol=1e-12)
        angle = np.deg2rad(prior['rotation_radius_degrees']); shift = prior['translation_radius_A']/g['field_A']
        grams, means = [], []
        for j in range(2):
            k, translations = perturbed_geometry(g['k'], g['q'], poses[j], angle, shift)
            gram = PlannedQuadratureObservationGram(k, g['ctf'], noise, order=40, preconditioner_rank=1024)
            wrapped = PhaseRotatedGram(gram, translations); grams.append(wrapped)
            means.append(wrapped.rotate(cell_forward(k, g['ctf'], pilot, 24, noise)))
        centers = [[0, 0, 0]] if target == 'center' else [[0, 0, .08], [0, 0, -.08]]
        signs = [1] if target == 'center' else [1, -1]
        fit = optimize_two_pose_modulus(grams, centers, signs, prior['width_fraction_field'], means, initial['weights'],
            prior['density_radius'], prior['pilot_norm'], prior['distance_budget'], maxiter=args.maxiter, callback=lambda r: print(r, flush=True))
        weights = fit.pop('weights'); upper = fit.pop('upper_weights'); witness = fit['witness']
        testing = float(2*norm.sf(witness['mean_distance_upper']/2))
        if testing <= 2*prior['alpha']:
            raise AssertionError('Testing condition failed')
        np.savez(path.with_suffix('.npz'), weights=weights, upper_weights=upper, pose0=poses[0], pose1=poses[1],
                 residual_amplitudes=witness['residual_amplitudes'], common_density_scale=witness['common_density_scale'],
                 indices=g['indices'], noise_std=noise)
        no_data = prior['density_radius']*fit['target_norm']
        np.testing.assert_allclose(no_data, prior['no_data_half_width'], rtol=1e-12)
        result.update(complete=True, numerical_failure=False, fit=fit, no_data_half_width=no_data,
                      minimum_sum_testing_error=testing, witness_array_sha256=sha(path.with_suffix('.npz')),
                      relative_half_width_lower=witness['fixed_length_half_width_lower']/no_data, seconds=time.perf_counter()-start,
                      previous_continuous_relative_half_width_lower=prior['relative_half_width_lower'],
                      projected_pose_step_relative_half_width_lower=pose_record['relative_half_width_lower'])
        save(); print('DONE', result['relative_half_width_lower'], fit['relative_modulus_gap'], flush=True)
    except Exception as exc:
        result.update(complete=False, error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
