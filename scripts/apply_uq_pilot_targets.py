#!/usr/bin/env python3
"""Conditional experimental application of every locked pilot-selected feature."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from fourier_splats.physics import fft_center
from fourier_splats.uq_data import particle_geometry, particle_observations, VoxelReference
from fourier_splats.uq_continuous import cell_forward, cell_target_coefficients
from fourier_splats.uq_noise_calibration import uncenter_fourier_weights
from fourier_splats.uq_group_noise import grouped_estimator_variance_upper, cell_density_distance
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_pilot_targets import read_target_lock, verify_pilot_scores, LOCK_SHA256
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'
MAPS = {'10028': '2660', '10049': '6487', '10076': '8434'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def completed(path):
    result = json.loads(path.read_text())
    if not result.get('complete') or result.get('error'):
        raise ValueError(f'Completed successful source required: {path}')
    return result


def run(args, locked, snapshot):
    dataset = locked['dataset']; out = BASE/args.output; out.mkdir(parents=True, exist_ok=True)
    original_path = BASE/'experimental-noise-grouped'/f'{dataset}.json'
    original = completed(original_path)
    selection = ROOT/f'research/uncertainty/confirmation/prediction-v1/{dataset}-selection.csv'
    if sha(selection) != original['source_selection_sha256']:
        raise ValueError('Fresh source selection changed')
    fresh_rows = list(csv.DictReader(selection.open()))
    old_rows = list(csv.DictReader((ROOT/f'research/uncertainty/splits/{dataset}.csv').open()))
    old_groups = {r['source_group'] for r in old_rows}
    indices = np.array(original['representative_indices'][original['calibration_groups']:])
    calibration_groups = np.array([fresh_rows[i]['source_group'] for i in indices])
    if len(set(calibration_groups)) != len(indices) or set(calibration_groups)&old_groups:
        raise ValueError('Require distinct fresh calibration groups disjoint from all old data')
    g = particle_geometry(ROOT, dataset, 'inference_half0', radius=12, count=128, seed=609315)
    cp = BASE/'representation'/dataset/'real_particles-spacing-2.0.npz'
    if sha(cp) != locked['pilot_checkpoint_sha256']:
        raise ValueError('Pilot changed')
    op, pilot, _, _ = model(g, np.load(cp), 24, noise=1.)
    pilot = op.expand(pilot); verify_pilot_scores(locked, pilot)
    amplitude = original['pilot_amplitude_raw_units']
    image_path = ROOT/f'data/uncertainty/confirmation/prediction-v1/{dataset}/images.npy'
    images = np.load(image_path, mmap_mode='r'); q = g['q'][0].astype(int); c = images.shape[-1]//2
    ft = fft_center(np.asarray(images[indices], float))[:, q[:, 1]+c, q[:, 0]+c]
    calibration = np.concatenate([ft.real, ft.imag], axis=1)/amplitude
    raw_observation = particle_observations(ROOT, dataset, g)
    y = np.concatenate([raw_observation.real, raw_observation.imag], axis=1).ravel()/amplitude
    phases = g['manifest']['data_sign']*np.exp(-2j*np.pi*np.einsum('nqi,ni->nq', g['q'], g['translations']))
    rp = ROOT/f'data/uncertainty/references/emd_{MAPS[dataset]}.map'
    reference = VoxelReference.from_mrc(rp, box=64).volume.ravel(); reference /= np.linalg.norm(reference)
    reference *= original['approx_reference_amplitude_raw_units']/amplitude
    distance, distance_pad = cell_density_distance(pilot, 24, reference, 64)
    beta = .005/12; delta = 1e-6/12; alpha_noise = (.045-1e-6)/12
    requested = args.features.split(',') if args.features else [f['name'] for f in locked['features']]
    if set(requested)-{f['name'] for f in locked['features']}:
        raise ValueError('Unknown target')
    for feature in locked['features']:
        name = feature['name']
        if name not in requested:
            continue
        path = out/f'{dataset}-{name}.json'
        if path.exists():
            if args.resume and completed(path)['target_lock_sha256'] == LOCK_SHA256:
                print('REUSE', path.name, flush=True); continue
            raise RuntimeError('Preserve previous outcome')
        start = time.perf_counter()
        fixed_path = BASE/args.fixed_source/f'{dataset}-{name}.json'; fixed = completed(fixed_path)
        if fixed['target_lock_sha256'] != LOCK_SHA256:
            raise ValueError('Fit does not use locked targets')
        row = next(r for r in fixed['targets'] if r['target'] == name); fit = row['fit']
        width = locked['width_fraction_field']; centers = [feature['center_fraction_field']]
        np.testing.assert_array_equal(row['centers_fraction_field'], centers)
        np.testing.assert_allclose(row['width_fraction_field'], width, rtol=1e-14)
        if row['signs'] != [1] or fixed['density_radius'] != 2. or fixed['pilot_norm_bound'] != 1.:
            raise ValueError('Source class changed')
        wp = fixed_path.with_name(f'{dataset}-{name}-{width}-weights.npz')
        if sha(wp) != row['source_weights_sha256']:
            raise ValueError('Source weights changed')
        saved = np.load(wp); np.testing.assert_array_equal(saved['indices'], g['indices'])
        w = saved['weights']; noise = float(saved['noise_std'])
        result = {'complete': False, 'dataset': dataset, 'target': name, 'width_A': 20.,
                  'centers_fraction_field': centers, 'target_lock_sha256': LOCK_SHA256, 'config': vars(args),
                  'scope': 'Exploratory conditional experimental sensitivity, not established density coverage',
                  'source_snapshot': snapshot, 'fixed_fit_sha256': sha(fixed_path), 'source_weights_sha256': sha(wp),
                  'source_experimental_amplitude_sha256': sha(original_path), 'reference_sha256': sha(rp),
                  'calibration_selection_sha256': sha(selection), 'calibration_indices': indices.tolist(),
                  'calibration_groups': calibration_groups.tolist(), 'old_calibration_group_overlap': 0,
                  'inference_indices': g['indices'].tolist(), 'inference_groups': g['groups'].tolist(),
                  'density_radius': 2., 'pilot_norm_bound': 1., 'pilot_amplitude_raw_units': amplitude,
                  'alpha_noise': alpha_noise, 'beta_calibration': beta, 'delta_spectral': delta,
                  'alpha_total_per_feature': .05/12, 'family_size': 12,
                  'reference_pilot_L2_distance': distance, 'reference_distance_squared_roundoff_pad': distance_pad,
                  'reference_inside_density_ball_numerical': bool(distance <= 2.),
                  'unverified_conditions': original['unverified_conditions']+[
                      'Inference exposures may contain arbitrary jointly Gaussian dependent particles; exposures remain independent and share one marginal raw noise covariance.',
                      'Inference assumes a shared density. In particular, the heterogeneous EMPIAR-10076 population does not establish this condition.',
                      'Calibration pool was examined in previous exploratory analyses; this is not new confirmation.'], 'records': []}
        def save():
            path.write_text(json.dumps(result, indent=2)+'\n')
        save()
        try:
            raw_weights = uncenter_fourier_weights(w.reshape(len(g['indices']), -1), phases, noise)
            calibrated = grouped_estimator_variance_upper(calibration, raw_weights, g['groups'], beta)
            sd = calibrated['noise_sd_upper']; result['group_noise_calibration'] = calibrated
            nominal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
            pilot_target = feature['pilot_expected_feature']
            observed = float(pilot_target+w@(y/noise-nominal))
            ref_target = float(cell_target_coefficients(64, centers, [1], width)@reference)
            no_data = 2*fit['target_norm']
            cases = [('fixed_pose', 0., 0., fit['bias'], fixed_path, None)]
            for angle in [0, 1, 2]:
                ap = BASE/args.pose_source/f'{dataset}-{name}-{angle}.json'; audit = completed(ap)
                ep = BASE/args.enclosure_source/ap.name; enclosure = completed(ep)
                if audit['source_fit_sha256'] != sha(fixed_path) or audit['source_weights_sha256'] != sha(wp):
                    raise ValueError('Pose source does not match fixed weights')
                if enclosure['source_audit_sha256'] != sha(ap) or enclosure['source_weights_sha256'] != sha(wp):
                    raise ValueError('Enclosure source changed')
                if audit['density_radius'] != 2. or audit['pilot_norm_bound'] != 1. or not enclosure['joint_ball']:
                    raise ValueError('Pose source has a different class')
                np.testing.assert_allclose(audit['alpha_numerical'], delta, rtol=1e-14)
                np.testing.assert_allclose([audit['config']['angle'], audit['config']['shift_A']], [angle, .5])
                cases.append(('shift_only' if angle == 0 else f'pose_{angle}deg', angle, .5, enclosure['bias_upper'], ep, ap))
            for label, angle, shift, bias, source, spectral_source in cases:
                half = bias_aware_half_width_stable(sd, bias, alpha_noise); fallback = bool(half >= no_data)
                center = pilot_target if fallback else observed; selected_half = min(half, no_data)
                result['records'].append({'pose_class': label, 'rotation_radius_degrees': angle, 'translation_radius_A': shift,
                    'bias_upper': bias, 'raw_observed_center': observed, 'noise_sd_upper': sd,
                    'raw_half_width': half, 'interval_center': center, 'interval_half_width': selected_half,
                    'relative_half_width': selected_half/no_data, 'uses_no_data': fallback,
                    'excludes_zero': bool(abs(center) > selected_half), 'approximate_reference_target': ref_target,
                    'approximate_reference_inside_interval': bool(abs(ref_target-center) <= selected_half),
                    'bias_source': str(source.relative_to(ROOT)), 'bias_source_sha256': sha(source),
                    'spectral_source': str(spectral_source.relative_to(ROOT)) if spectral_source else None,
                    'spectral_source_sha256': sha(spectral_source) if spectral_source else None})
            result.update(complete=True, seconds=time.perf_counter()-start); save()
            print('DONE', dataset, name, [(r['pose_class'], r['relative_half_width'], r['excludes_zero']) for r in result['records']], flush=True)
        except Exception as exc:
            result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


def main():
    p = argparse.ArgumentParser(); p.add_argument('--datasets', default='10028,10049,10076'); p.add_argument('--features')
    p.add_argument('--fixed-source', default='pilot-selected-fixed-v1'); p.add_argument('--pose-source', default='pilot-selected-pose-v1')
    p.add_argument('--enclosure-source', default='pilot-selected-enclosing-v1'); p.add_argument('--resume', action='store_true')
    p.add_argument('--output', default='pilot-selected-experimental-v1'); args = p.parse_args()
    lock = read_target_lock(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
    if set(args.datasets.split(','))-{r['dataset'] for r in lock['datasets']}:
        raise ValueError('Unknown dataset')
    snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py',
        'research/uncertainty/pilot-selected-targets-v1/EXPERIMENTAL-APPLICATION.md'])
    for row in lock['datasets']:
        if row['dataset'] in args.datasets.split(','):
            run(args, row, snapshot)


if __name__ == '__main__':
    main()
