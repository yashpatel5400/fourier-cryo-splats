#!/usr/bin/env python3
"""Fixed-pose Fourier VI and exact Gaussian posterior on three stack geometries.

Known simulated noise, independently specified prior scales, and explicit
continuous-map/model-matched/pose-perturbed controls. Not an external-code run
or empirical coverage estimate for raw cryo-EM density.
"""
import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry, VoxelReference
from fourier_splats.uq_continuous import cell_target_coefficients, cell_forward
from fourier_splats.uq_continuous_pose import pose_cell_forward
from fourier_splats.uq_fourier_variational import HermitianTrilinearOperator, diagonal_variational_variances
from fourier_splats.uq_baselines import gaussian_reference_operator
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'
MAPS = {'10028': '2660', '10049': '6487', '10076': '8434'}


def run(args, dataset, snapshot):
    out = BASE/args.output; out.mkdir(exist_ok=True); path = out/f'{dataset}.json'
    if path.exists(): raise RuntimeError('Preserve previous baseline outcomes')
    start = time.perf_counter()
    record = {'stage': 'Equation-14 Gaussian VI reimplementation with CTFs and declared priors; fixed experimental geometry and known simulated noise',
              'dataset': dataset, 'config': vars(args), 'source_snapshot': snapshot,
              'complete': False, 'records': [], 'source': 'https://proceedings.mlr.press/v115/ullrich20a.html'}
    def save(): path.write_text(json.dumps(record, indent=2)+'\n')
    save()
    try:
        g = particle_geometry(ROOT, dataset, 'inference_half0', radius=args.radius, count=args.particles, seed=609841)
        noise = float(np.load(BASE/'continuous-quadrature-optimized'/f'{dataset}-center-0.07-weights.npz')['noise_std'])
        ck = np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
        old_op, pilot, _, _ = model(g, ck, 24, noise=noise); pilot = old_op.expand(pilot)
        ref = VoxelReference.from_mrc(ROOT/f'data/uncertainty/references/emd_{MAPS[dataset]}.map', box=64).volume.ravel()
        ref /= np.linalg.norm(ref)
        op = HermitianTrilinearOperator(g['k'], g['ctf'], noise, box=args.box)
        pilot_x = op.project_cells(pilot, 24); reference_x = op.project_cells(ref, 64)
        pilot_signal = op.matrix@pilot_x; matched_signal = op.matrix@reference_x
        exact_signal = cell_forward(g['k'], g['ctf'], ref, 64, noise)
        signals = [('matched_fourier_model', matched_signal), ('continuous_reference_nominal', exact_signal)]
        rng = np.random.default_rng(609851+int(dataset))
        for degrees in [.5, 1., 2.]:
            for scenario in ['coherent_x', 'random_boundary']:
                u = np.zeros((len(g['k']), 5))
                if scenario == 'coherent_x': u[:, 0] = 1.
                else:
                    u = rng.normal(size=u.shape); u /= np.linalg.norm(u, axis=1)[:, None]
                signal = pose_cell_forward(g['k'], g['q'], g['ctf'], ref, 64, noise, u,
                                           np.deg2rad(degrees), .5/g['field_A'])
                signals.append((f'continuous_reference_{degrees:g}deg_{scenario}', signal))
        record.update(noise_std_supplied=noise, field_A=g['field_A'], particles=len(g['k']),
            independent_fourier_pairs_per_particle=g['k'].shape[1], parameters=op.shape[1],
            sparse_nonzeros=int(op.matrix.nnz), sparse_matrix_bytes=int(op.matrix.data.nbytes+op.matrix.indices.nbytes+op.matrix.indptr.nbytes),
            geometry_indices_sha256=hashlib.sha256(g['indices'].tobytes()).hexdigest(),
            continuous_reference_relative_forward_error=float(np.linalg.norm(matched_signal-exact_signal)/np.linalg.norm(exact_signal)),
            prior_center_norm=float(np.linalg.norm(pilot_x)), projected_reference_norm=float(np.linalg.norm(reference_x)),
            setup_seconds=time.perf_counter()-start)
        save()
        for target in ['center', 'contrast']:
            offset = args.width_A/g['field_A']
            centers = [[0, 0, 0]] if target == 'center' else [[0, 0, offset], [0, 0, -offset]]
            signs = [1] if target == 'center' else [1, -1]
            ell, full_norm2 = op.target(centers, signs, offset)
            full_truth = float(cell_target_coefficients(64, centers, signs, offset)@ref)
            matched_truth = float(ell@reference_x); prior_target = float(ell@pilot_x)
            for prior_norm in map(float, args.prior_norms.split(',')):
                tick = time.perf_counter(); tau = prior_norm/np.sqrt(op.shape[1])
                fit = gaussian_reference_operator(op.matrix, ell, tau, gram_diagonal=op.gram_diagonal, rtol=1e-10)
                w = fit.pop('weights'); residual = fit.pop('density_bias_direction')
                vi_var = float(ell@(diagonal_variational_variances(op.gram_diagonal, tau)*ell))
                vi_half = float(norm.ppf(.975)*np.sqrt(vi_var))
                marginal_var = float(w@w+tau**2*(residual@residual))
                post_var = (fit['posterior_half_width']/norm.ppf(.975))**2
                identity_error = abs(marginal_var-post_var)/max(post_var, 1e-300)
                row = {'target': target, 'sigma_A': args.width_A, 'contrast_centers_A': ([-args.width_A, args.width_A] if target == 'contrast' else [0]),
                    'prior_expected_squared_deviation_norm': prior_norm**2, 'prior_coordinate_sd': tau,
                    'full_posterior': fit, 'diagonal_vi_half_width': vi_half,
                    'vi_over_full_posterior_width': vi_half/fit['posterior_half_width'],
                    'continuous_target_norm2': full_norm2, 'model_target_norm2': float(ell@ell),
                    'prior_predictive_variance_identity_relative_error': identity_error,
                    'full_prior_predictive_coverage': float(2*norm.cdf(fit['posterior_half_width']/np.sqrt(marginal_var))-1),
                    'diagonal_vi_prior_predictive_coverage': float(2*norm.cdf(vi_half/np.sqrt(marginal_var))-1),
                    'checks': []}
                sd = float(np.linalg.norm(w))
                for scenario, signal in signals:
                    truth = matched_truth if scenario == 'matched_fourier_model' else full_truth
                    mean = float(prior_target+w@(signal-pilot_signal)); bias = mean-truth
                    for method, half in [('diagonal_fourier_vi', vi_half), ('full_fourier_gaussian', fit['posterior_half_width'])]:
                        coverage = float(norm.cdf((half-bias)/sd)-norm.cdf((-half-bias)/sd)) if sd else float(abs(bias) <= half)
                        power = float(norm.cdf((np.sign(truth)*mean-half)/sd)) if sd else float(np.sign(truth)*mean > half)
                        row['checks'].append({'scenario': scenario, 'method': method, 'true_target': truth,
                            'expected_center': mean, 'bias': bias, 'half_width': half,
                            'analytic_fixed_signal_coverage': coverage, 'correct_sign_probability': power})
                row['seconds'] = time.perf_counter()-tick
                row['numerical_failure'] = bool(identity_error > 1e-7)
                record['records'].append(row); save()
                np.savez(out/f'{dataset}-{target}-prior{prior_norm:g}.npz', weights=w, indices=g['indices'],
                         noise_std=noise, target_coefficients=ell, prior_center_coefficients=pilot_x)
                print(dataset, target, prior_norm, 'vi/full', row['vi_over_full_posterior_width'],
                      'VI prior coverage', row['diagonal_vi_prior_predictive_coverage'],
                      'identity error', identity_error, flush=True)
                if row['numerical_failure']: raise AssertionError('Posterior identity check failed; case retained')
        record.update(complete=True, seconds=time.perf_counter()-start,
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024))
        save()
    except Exception as exc:
        record.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--datasets', default='10028,10049,10076')
    p.add_argument('--particles', type=int, default=1024); p.add_argument('--radius', type=float, default=12)
    p.add_argument('--box', type=int, default=33); p.add_argument('--width-A', type=float, default=10)
    p.add_argument('--prior-norms', default='.5,1,2'); p.add_argument('--output', default='fourier-variational-baseline')
    args = p.parse_args(); snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py'])
    for dataset in args.datasets.split(','): run(args, dataset, snapshot)
