#!/usr/bin/env python3
"""Declared local Gaussian pose-marginalized baseline for the locked target family."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import resource
import sys
import subprocess
from scipy.spatial.transform import Rotation
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry, VoxelReference
from fourier_splats.uq_continuous import cell_target_coefficients, cell_forward
from fourier_splats.uq_continuous_pose import pose_cell_forward
from fourier_splats.uq_fourier_variational import HermitianTrilinearOperator
from fourier_splats.uq_fourier_pose_baseline import trilinear_pose_jacobian, GaussianNuisanceWhitening
from fourier_splats.uq_baselines import gaussian_reference_operator
from fourier_splats.uq_pilot_targets import read_target_lock, verify_pilot_scores, LOCK_SHA256
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'
MAPS = {'10028': '2660', '10049': '6487', '10076': '8434'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(locked, snapshot, args):
    dataset = locked['dataset']; out = BASE/args.output; out.mkdir(exist_ok=True)
    path = out/f'{dataset}.json'
    if path.exists():
        raise RuntimeError('Preserve previous baseline outcomes')
    start = time.perf_counter(); alpha = .05/12; z = norm.ppf(1-alpha/2)
    scales = [.1, 1.]
    if not scales or any(not np.isfinite(s) or s <= 0 for s in scales):
        raise ValueError('Finite positive prior scales required')
    result = {'complete': False, 'dataset': dataset, 'target_lock_sha256': LOCK_SHA256,
        'source_snapshot': snapshot, 'alpha_per_feature': alpha, 'family_size_per_prior_scale': 12,
        'scope': 'Local Gaussian pose marginalization around a fixed pilot; supplied priors and known noise, not end-to-end experimental calibration',
        'source': 'https://proceedings.mlr.press/v115/ullrich20a.html',
        'config': {'particles': 128, 'radius': 12, 'geometry_seed': 609315, 'box': 33,
                   'prior_coordinate_sds': scales, 'pose_models': ['fixed_pose', 'local_gaussian_pose'], 'reference_pose_seed': 610281+int(dataset),
                   'prior_scale_source': 'Previously inspected broader priors, retained together', 'pose_coordinate_sd_degrees': 1., 'shift_coordinate_sd_A': .5}, 'records': []}
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    save()
    try:
        g = particle_geometry(ROOT, dataset, 'inference_half0', radius=12, count=128, seed=609315)
        noise_path = BASE/'continuous-quadrature-optimized'/f'{dataset}-center-0.07-weights.npz'
        noise = float(np.load(noise_path)['noise_std'])
        cp = BASE/f'representation/{dataset}/real_particles-spacing-2.0.npz'
        if sha(cp) != locked['pilot_checkpoint_sha256']:
            raise ValueError('Pilot differs from target lock')
        old_op, pilot, _, _ = model(g, np.load(cp), 24, noise=noise)
        pilot = old_op.expand(pilot); verify_pilot_scores(locked, pilot)
        rp = ROOT/f'data/uncertainty/references/emd_{MAPS[dataset]}.map'
        ref = VoxelReference.from_mrc(rp, box=64).volume.ravel(); ref /= np.linalg.norm(ref)
        op = HermitianTrilinearOperator(g['k'], g['ctf'], noise, box=33)
        pilot_x = op.project_cells(pilot, 24); ref_x = op.project_cells(ref, 64)
        pilot_signal = op.matrix@pilot_x; exact = cell_forward(g['k'], g['ctf'], ref, 64, noise)
        matched = op.matrix@ref_x
        jacobian, diagnostics = trilinear_pose_jacobian(op, g['k'], g['q'], g['ctf'], noise, pilot_x)
        scaled_jacobian = jacobian*np.array([np.deg2rad(1.)]*3+[.5/g['field_A']]*2)
        whitener = GaussianNuisanceWhitening(scaled_jacobian)
        model_operators = [('fixed_pose', op.matrix, None),
                           ('local_gaussian_pose', whitener.whiten_operator(op.matrix), whitener)]
        jp = out/f'{dataset}-jacobian.npz'
        np.savez(jp, scaled_jacobian=scaled_jacobian, indices=g['indices'])
        result.update(jacobian_sha256=sha(jp), jacobian_diagnostics=diagnostics, linearization_checks=[])
        old_path = BASE/f'pilot-selected-fourier-baselines-weak/{dataset}.json'
        old = json.loads(old_path.read_text())
        if not old.get('complete') or old.get('error'):
            raise ValueError('Completed original weak-prior baseline required')
        result['original_fixed_baseline_sha256'] = sha(old_path)
        signals = [('matched_fourier_model', matched), ('nominal', exact)]
        for angle in [0., 1., 2.]:
            rng = np.random.default_rng(610281+int(dataset))
            for scenario in ['coherent_x', 'random_boundary']:
                u = np.zeros((128, 5))
                if scenario == 'coherent_x':
                    u[:, 0] = 1.
                else:
                    u = rng.normal(size=u.shape); u /= np.linalg.norm(u, axis=1)[:, None]
                signal = pose_cell_forward(g['k'], g['q'], g['ctf'], ref, 64, noise, u,
                                           np.deg2rad(angle), .5/g['field_A'])
                signals.append((f'{angle:g}deg_{scenario}', signal))
                rotated = g['k']@Rotation.from_rotvec(np.deg2rad(angle)*u[:, :3]).as_matrix()
                changed = (HermitianTrilinearOperator(rotated, g['ctf'], noise, box=33).matrix@pilot_x).reshape(128, -1)
                nq = g['q'].shape[1]
                phase = np.exp(-2j*np.pi*.5/g['field_A']*np.einsum('nqi,ni->nq', g['q'], u[:, 3:]))
                complex_value = (changed[:, :nq]+1j*changed[:, nq:])*phase
                nonlinear_pilot = np.concatenate([complex_value.real, complex_value.imag], axis=1).ravel()
                physical_pose = u*np.array([np.deg2rad(angle)]*3+[.5/g['field_A']]*2)
                linear_pilot = pilot_signal+np.einsum('nmp,np->nm', jacobian, physical_pose).ravel()
                delta = nonlinear_pilot-pilot_signal; error = nonlinear_pilot-linear_pilot
                result['linearization_checks'].append({'scenario': f'{angle:g}deg_{scenario}',
                    'absolute_error': float(np.linalg.norm(error)),
                    'relative_to_signal': float(np.linalg.norm(error)/np.linalg.norm(nonlinear_pilot)),
                    'relative_to_pose_change': float(np.linalg.norm(error)/max(np.linalg.norm(delta),1e-300))})
        result.update(noise_std_supplied=noise, noise_source_sha256=sha(noise_path),
            field_A=g['field_A'], band_endpoint_A=g['field_A']/12,
            geometry_indices_sha256=hashlib.sha256(g['indices'].tobytes()).hexdigest(),
            reference_sha256=sha(rp), pilot_checkpoint_sha256=sha(cp),
            forward_interpolation_relative_error=float(np.linalg.norm(matched-exact)/np.linalg.norm(exact)),
            prior_center_norm=float(np.linalg.norm(pilot_x)), setup_seconds=time.perf_counter()-start)
        save()
        for feature in locked['features']:
            name = feature['name']; centers = [feature['center_fraction_field']]
            width = locked['width_fraction_field']; ell, norm2 = op.target(centers, [1], width)
            truth = float(cell_target_coefficients(64, centers, [1], width)@ref)
            matched_truth = float(ell@ref_x); prior_target = float(ell@pilot_x)
            no_data = 2*np.sqrt(norm2)
            for tau in scales:
                for pose_model, design, transform in model_operators:
                    tick = time.perf_counter()
                    fit = gaussian_reference_operator(design, ell, tau, alpha=alpha,
                        gram_diagonal=op.gram_diagonal, rtol=1e-10, maxiter=2000)
                    fitted_weights = fit.pop('weights'); residual = fit.pop('density_bias_direction')
                    w = fitted_weights if transform is None else transform.apply(fitted_weights)
                    full_half = fit['posterior_half_width']; sampling_sd = float(np.linalg.norm(w))
                    fit['measurement_noise_half_width'] = float(z*sampling_sd)
                    pose_variance = 0. if transform is None else float(np.sum(np.einsum('nmp,nm->np', scaled_jacobian, w.reshape(128,-1))**2))
                    prior_variance = float(w@w+pose_variance+tau**2*(residual@residual))
                    posterior_variance = (full_half/z)**2
                    error = abs(prior_variance-posterior_variance)/max(posterior_variance, 1e-300)
                    previous = next(r for r in old['records'] if r['target']==name and abs(r['prior_coordinate_sd']-tau)<1e-12)
                    old_half = previous['full_posterior']['posterior_half_width']
                    fixed_repeat_error = abs(full_half-old_half)/old_half if transform is None else None
                    monotonicity_failure = bool(full_half < old_half*(1-2e-8))
                    fixed_repeat_failure = bool(transform is None and fixed_repeat_error > 2e-8)
                    wp = out/f'{dataset}-{name}-prior{tau:g}-{pose_model}.npz'
                    np.savez(wp, weights=w, indices=g['indices'], noise_std=noise,
                        target_coefficients=ell, prior_center_coefficients=pilot_x)
                    row = {'target': name, 'centers_fraction_field': centers, 'width_A': 20.,
                        'prior_deviation_norm': tau*np.sqrt(op.shape[1]), 'prior_coordinate_sd': tau,
                        'pose_model': pose_model, 'full_posterior': fit,
                        'continuous_target_norm': float(np.sqrt(norm2)), 'model_target_norm': float(np.linalg.norm(ell)),
                        'continuous_pilot_feature': feature['pilot_expected_feature'], 'projected_prior_feature': prior_target,
                        'width_over_original_fixed': full_half/old_half,
                        'prior_variance_identity_relative_error': error,
                        'fixed_repeat_relative_error': fixed_repeat_error,
                        'variance_monotonicity_failure': monotonicity_failure,
                        'numerical_failure': bool(error>1e-7 or monotonicity_failure or fixed_repeat_failure),
                        'measurement_variance': float(w@w), 'local_pose_variance': pose_variance,
                        'prior_bias_variance': float(tau**2*(residual@residual)),
                        'full_prior_predictive_coverage': float(2*norm.cdf(full_half/np.sqrt(prior_variance))-1),
                        'weights_sha256': sha(wp), 'checks': []}
                    for scenario, signal in signals:
                        actual = matched_truth if scenario == 'matched_fourier_model' else truth
                        mean = float(prior_target+w@(signal-pilot_signal)); bias = mean-actual
                        coverage = norm.cdf((full_half-bias)/sampling_sd)-norm.cdf((-full_half-bias)/sampling_sd)
                        power = norm.cdf((np.sign(actual)*mean-full_half)/sampling_sd)
                        row['checks'].append({'scenario': scenario,
                            'true_target': actual, 'expected_center': mean, 'bias': bias,
                            'measurement_noise_sd': sampling_sd, 'half_width': full_half, 'relative_half_width': full_half/no_data,
                            'analytic_fixed_signal_coverage': float(coverage), 'correct_sign_probability': float(power)})
                    row['seconds'] = time.perf_counter()-tick
                    result['records'].append(row); save()
                    if row['numerical_failure']:
                        raise ArithmeticError('Posterior identity, variance ordering or fixed-control repeat failed; full outcome preserved')
                    print('DONE',dataset,name,tau,pose_model,'pose/fixed width',row['width_over_original_fixed'],flush=True)
        result.update(complete=True, seconds=time.perf_counter()-start, peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)); save()
    except Exception as exc:
        result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


def main():
    p = argparse.ArgumentParser(); p.add_argument('--datasets', default='10028,10049,10076')
    p.add_argument('--output', default='pilot-selected-fourier-pose-v1'); args = p.parse_args()
    lock = read_target_lock(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
    if set(args.datasets.split(','))-{d['dataset'] for d in lock['datasets']}:
        raise ValueError('Unknown dataset')
    for file in [Path(__file__), *sorted((ROOT/'src/fourier_splats').glob('*.py')), ROOT/'research/uncertainty/FOURIER-POSE-BASELINE-PROTOCOL.md']:
        if subprocess.check_output(['git','show',f'HEAD:{file.relative_to(ROOT)}'],cwd=ROOT) != file.read_bytes():
            raise ValueError('Commit baseline source before fitting')
    snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py',
        'research/uncertainty/pilot-selected-targets-v1/BASELINES.md',
        'research/uncertainty/FOURIER-POSE-BASELINE-PROTOCOL.md', 'tests/test_uq_fourier_pose_baseline.py',
        'research/uncertainty/pilot-selected-targets-v1/locked-targets.json'])
    errors = []
    for row in lock['datasets']:
        if row['dataset'] in args.datasets.split(','):
            try:
                run(row, snapshot, args)
            except Exception as exc:
                errors.append((row['dataset'],repr(exc))); print('RETAINED_FAILURE',errors[-1],flush=True)
    if errors:
        raise RuntimeError(f'All requested datasets attempted; failures: {errors}')


if __name__ == '__main__':
    main()
