#!/usr/bin/env python3
"""Feasible pose refinement of one projected constant-cell ambiguity pair."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous import cell_adjoint, cell_target_coefficients
from fourier_splats.uq_continuous_pose import perturbed_geometry
from fourier_splats.uq_cell_pose_pair import CellPairDistance, optimize_cell_pair_poses, refine_cell_density
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', default='two-pose-modulus-v3/10049-center-random_boundary-2.json')
    p.add_argument('--max-evaluations', type=int, default=80)
    p.add_argument('--box', type=int, choices=[24, 48, 72], default=24)
    p.add_argument('--output', default='pose-optimized-ambiguity-pilot'); args = p.parse_args()
    source = BASE/args.source; prior = json.loads(source.read_text()); wp = source.with_suffix('.npz')
    if not prior.get('complete') or prior.get('numerical_failure') or prior.get('error'):
        raise ValueError('Completed feasible source witness required')
    if prior['pose_set'] != 'joint_scaled_five_dimensional_ball' or prior['density_radius'] != 2. or abs(prior['pilot_norm']-1) > 1e-12:
        raise ValueError('This protocol requires the original joint pose ball, B=2 and unit pilot')
    if sha(wp) != prior['witness_array_sha256']:
        raise ValueError('Source witness changed')
    out = BASE/args.output; out.mkdir(parents=True, exist_ok=True); path = out/source.name
    if path.exists():
        raise RuntimeError('Preserve previous outcome')
    start = time.perf_counter()
    result = {'complete': False, 'scope': 'Constructive conditional constant-cell ambiguity pair with locally optimized poses',
              'config': vars(args), 'source_record_sha256': sha(source), 'source_arrays_sha256': sha(wp),
              'source_snapshot': source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py', 'research/uncertainty/POSE-OPTIMIZED-AMBIGUITY-PROTOCOL.md']),
              'dataset': prior['dataset'], 'target': prior['target'], 'alpha': prior['alpha'],
              'density_radius': prior['density_radius'], 'pilot_norm': prior['pilot_norm'],
              'pose_set': prior['pose_set'], 'rotation_radius_degrees': prior['rotation_radius_degrees'],
              'translation_radius_A': prior['translation_radius_A'], 'width_fraction_field': prior['width_fraction_field'],
              'no_data_half_width': prior['no_data_half_width'], 'source_relative_half_width_lower': prior['relative_half_width_lower']}
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    save()
    try:
        dataset = prior['dataset']; target = prior['target']; box = args.box; B = prior['density_radius']
        g = particle_geometry(ROOT, dataset, 'inference_half0', radius=5, count=128, seed=prior['source_geometry_config']['seed'])
        saved = np.load(wp); np.testing.assert_array_equal(saved['indices'], g['indices']); noise = float(saved['noise_std'])
        checkpoint = BASE/'representation'/dataset/'real_particles-spacing-2.0.npz'
        if sha(checkpoint) != prior['pilot_checkpoint_sha256']:
            raise ValueError('Pilot changed')
        cellop, pilot, _, check_noise = model(g, np.load(checkpoint), 24); pilot = cellop.expand(pilot)
        pilot = refine_cell_density(pilot, 24, box)
        np.testing.assert_allclose(noise, check_noise); np.testing.assert_allclose(np.linalg.norm(pilot), 1., rtol=1e-12)
        centers = [[0, 0, 0]] if target == 'center' else [[0, 0, .08], [0, 0, -.08]]
        signs = [1] if target == 'center' else [1, -1]
        ell = cell_target_coefficients(box, centers, signs, prior['width_fraction_field'])
        n, nq = g['ctf'].shape; w = saved['weights'].reshape(n, 2*nq)
        angle = np.deg2rad(prior['rotation_radius_degrees']); shift = prior['translation_radius_A']/g['field_A']
        poses = np.stack([saved['pose0'], saved['pose1']]); densities = []; inward = []
        for j in range(2):
            k, phase = perturbed_geometry(g['k'], g['q'], poses[j], angle, shift)
            rotated_weights = (w[:, :nq]+1j*w[:, nq:])*np.exp(2j*np.pi*phase)
            weights = np.concatenate([rotated_weights.real, rotated_weights.imag], axis=1).ravel()
            h = ell-cell_adjoint(k, g['ctf'], weights, box, noise)
            delta = (-1 if j == 0 else 1)*float(saved['residual_amplitudes'][j])*h
            factor = min(1., B*(1-1e-10)/max(np.linalg.norm(delta), np.finfo(float).tiny))
            inward.append(float(factor)); densities.append(pilot+factor*delta)
        densities = np.stack(densities)
        raw_gap = float(ell@(densities[1]-densities[0]))
        gap_pad = 128*np.finfo(float).eps*max(1., np.sum(abs(ell*(densities[1]-densities[0]))))
        gap = max(0., abs(raw_gap)-gap_pad); tau = prior['distance_budget']
        problem = CellPairDistance(g['k'], g['q'], g['ctf'], densities, box, noise, angle, shift)
        initial_distance, initial_check = problem.independent_distance(poses)
        initial_scale = min(1., tau/initial_distance)
        result.update(box=box, noise_std=noise, field_A=g['field_A'], inward_projection_factors=inward,
                      projected_class_distances=np.linalg.norm(densities-pilot, axis=1).tolist(),
                      projected_unscaled_feature_gap=gap, target_gap_downward_pad=float(gap_pad),
                      projected_initial_distance=initial_distance, projected_initial_check=initial_check,
                      projected_initial_relative_half_width_lower=initial_scale*gap/2/prior['no_data_half_width'])
        save()
        fit = optimize_cell_pair_poses(problem, poses, max_evaluations=args.max_evaluations, callback=lambda r: print(r, flush=True))
        optimized = fit.pop('poses')/(1+1e-12)
        distance, final_check = problem.independent_distance(optimized)
        scale = min(1., tau/distance); final = scale*densities
        class_norms = np.linalg.norm(final-pilot, axis=1); pose_norm = float(np.max(np.linalg.norm(optimized, axis=-1)))
        if class_norms.max() > B or pose_norm > 1.:
            raise AssertionError('Final pair violates class')
        testing = float(2*norm.sf(scale*distance/2))
        if testing <= 2*prior['alpha']:
            raise AssertionError('Testing separation too large')
        np.savez(path.with_suffix('.npz'), densities=final, unscaled_densities=densities, pilot=pilot, poses=optimized,
                 initial_poses=poses, target_coefficients=ell, indices=g['indices'], noise_std=noise, common_scale=scale)
        result.update(complete=True, fit=fit, final_distance_check=final_check, mean_distance_upper=scale*distance,
                      unscaled_mean_distance_upper=distance, common_density_scale=scale,
                      final_class_distances=class_norms.tolist(), maximum_pose_norm=pose_norm,
                      final_negative_density_energy_fractions=(np.sum(np.minimum(final, 0.)**2, axis=1)/np.sum(final**2, axis=1)).tolist(),
                      minimum_sum_testing_error=testing, fixed_length_half_width_lower=scale*gap/2,
                      relative_half_width_lower=scale*gap/2/prior['no_data_half_width'],
                      arrays_sha256=sha(path.with_suffix('.npz')), seconds=time.perf_counter()-start,
                      numerical_scope='Exact cell Fourier integrals in real arithmetic; NUFFT and floating magnitude guards are not interval validated')
        save(); print('DONE', result['relative_half_width_lower'], result['projected_initial_relative_half_width_lower'], flush=True)
    except Exception as exc:
        result.update(complete=False, error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
