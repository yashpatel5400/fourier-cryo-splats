#!/usr/bin/env python3
"""One bounded cubic weight-design case, declared before empirical fitting."""
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import finufft
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_cubic_optimization import CubicDesignObjective, optimize_cubic_weights
from fourier_splats.uq_cubic_design import DifferentiableCubicPoseFieldOperator
from fourier_splats.uq_cubic_pose import cell_moments_cubic, cubic_residual_cross
from fourier_splats.uq_sobolev_penalty import SobolevRemainderPenalty
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_continuous import cell_forward, cell_target_coefficients
from fourier_splats.uq_continuous_pose import pose_cell_forward
from fourier_splats.uq_data import particle_geometry, VoxelReference
from fourier_splats.uq_higher_remainder import higher_order_particle_remainders
from fourier_splats.uq_joint_bias import joint_density_pose_bias
from fourier_splats.uq_intervals import bias_aware_half_width_stable, reference_interval_summary
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
SOURCE = BASE/'cubic-pose-probe/10049-pilot_region_1-1.json'
PROTOCOL = ROOT/'research/uncertainty/CUBIC-WEIGHT-OPTIMIZATION-PROTOCOL.md'
CONFIG = {'dataset': '10049', 'target': 'pilot_region_1', 'particles': 128, 'radius': 12,
          'width_A': 20., 'angle_degrees': 1., 'shift_A': .5, 'density_radius': 2., 'pilot_norm': 1.,
          'density_order': 80, 'pose_order': 64, 'preconditioner_rank': 1024, 'threads': 2,
          'maxiter': 30, 'max_evaluations': 40, 'maxls': 10, 'ftol': 1e-5, 'gtol': 1e-6,
          'modes': 3, 'ritz_tolerance': 1e-3, 'ritz_subspace': 13, 'ritz_maxiter': 100, 'smoothing': .002,
          'optimization_seed': 650011, 'certificate_seed': 650001, 'probes': 4, 'power_iterations': 40,
          'alpha_total': .05/12, 'delta': 1e-6/12}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seed_inventory():
    matches = []; count = 0
    def scan(value, keys, name):
        if isinstance(value, dict):
            for key, child in value.items():
                if 'seed' in key and child == CONFIG['certificate_seed']:
                    matches.append({'path': name, 'keys': keys+[key]})
                scan(child, keys+[key], name)
        elif isinstance(value, list):
            for j, child in enumerate(value):
                scan(child, keys+[j], name)
    for path in sorted(BASE.rglob('*.json')):
        count += 1; scan(json.loads(path.read_text()), [], str(path.relative_to(ROOT)))
    if matches:
        raise ValueError(f'Certificate seed already recorded: {matches}')
    return {'json_files_examined': count, 'matches': matches, 'seed': CONFIG['certificate_seed'],
            'scope': 'All existing development JSON seed-valued fields before creating this run; not text substring matches.'}


def main():
    if 'tmp/nufft-openmp' not in finufft.__file__ or 'torch' in sys.modules:
        raise RuntimeError('Isolated threaded CPU-only FINUFFT required')
    out = BASE/'cubic-weight-probe'; path = out/'10049-pilot_region_1-1.json'
    if path.exists():
        raise RuntimeError('Preserve existing attempt, including failures')
    inventory = seed_inventory()
    old = json.loads(SOURCE.read_text()); source_audit = ROOT/old['source_audit']
    original = json.loads(source_audit.read_text()); source_fit = ROOT/original['config']['fit']
    fit = json.loads(source_fit.read_text()); cfg = fit['config']; row = fit['targets'][0]
    if not all(r.get('complete') and not r.get('error') for r in (old, original, fit)):
        raise ValueError('Completed immutable source cases required')
    if sha(source_audit) != old['source_audit_sha256'] or sha(source_fit) != old['source_fit_sha256']:
        raise ValueError('Source case changed')
    if cfg['particles'] != 128 or cfg['frequency_radius'] != 12 or len(fit['targets']) != 1:
        raise ValueError('Changed empirical cohort')
    if old['dataset'] != '10049' or old['target'] != 'pilot_region_1' or old['density_radius'] != 2. or old['pilot_norm_bound'] != 1.:
        raise ValueError('Changed empirical target/class')
    width = original['width_fraction_field']; centers = original['centers_fraction_field']; signs = original['signs']
    wp = source_fit.with_name(f'10049-pilot_region_1-{width}-weights.npz')
    if sha(wp) != old['source_weights_sha256']:
        raise ValueError('Starting weights changed')
    extra = ['scripts/audit_uq_grid_refinement.py', 'research/uncertainty/CUBIC-WEIGHT-OPTIMIZATION-PROTOCOL.md',
             'tests/test_uq_cubic_optimization.py']
    # No empirical run may silently execute unpublished changed implementation.
    for p in [Path(__file__), *sorted((ROOT/'src/fourier_splats').glob('*.py')), *[ROOT/f for f in extra]]:
        relative = str(p.relative_to(ROOT))
        if subprocess.check_output(['git', 'show', f'HEAD:{relative}'], cwd=ROOT) != p.read_bytes():
            raise ValueError(f'Commit declared source before empirical fitting: {relative}')
    out.mkdir(parents=True, exist_ok=True)
    result = {'complete': False, 'config': CONFIG, 'dataset': '10049', 'target': 'pilot_region_1',
        'scope': 'Single declared development fit with known noise and assumed joint pose balls; not fresh experimental calibration.',
        'source_case': str(SOURCE.relative_to(ROOT)), 'source_case_sha256': sha(SOURCE),
        'source_fit': str(source_fit.relative_to(ROOT)), 'source_fit_sha256': sha(source_fit),
        'source_weights_sha256': sha(wp), 'protocol_sha256': sha(PROTOCOL), 'seed_inventory': inventory,
        'source_snapshot': source_snapshot(ROOT, Path(__file__), extra),
        'checkpoints': [], 'reference_checks': [], 'pose_scaling': old['pose_scaling']}
    start = time.perf_counter()
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    def progress(record):
        print(record, flush=True)
        result['last_progress'] = record; save()
    def checkpoint(weights, record):
        destination = out/f'10049-pilot_region_1-1-eval-{record["evaluation"]:04d}.npz'
        if destination.exists():
            raise RuntimeError('Checkpoint collision')
        np.savez(destination, weights=weights)
        result['checkpoints'].append({'path': str(destination.relative_to(ROOT)), 'sha256': sha(destination), 'evaluation': record})
        save()
    save()
    try:
        g = particle_geometry(ROOT, '10049', 'inference_half0', radius=12, count=128, seed=cfg['seed'])
        saved = np.load(wp); np.testing.assert_array_equal(saved['indices'], g['indices'])
        np.testing.assert_allclose(width*g['field_A'], 20., rtol=1e-12)
        w = saved['weights']; noise = float(saved['noise_std']); angle = np.deg2rad(1.); shift = .5/g['field_A']
        gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=80, preconditioner_rank=1024)
        gram.nthreads = 2
        progress({'stage': 'nominal_gram_ready', 'seconds': time.perf_counter()-start})
        op = DifferentiableCubicPoseFieldOperator(g['k'], g['q'], g['ctf'], w, noise, angle, shift,
            order=64, nthreads=2, block_particles=16)
        scales = np.asarray(old['pose_scaling']['group_scales'], dtype='<f8')
        denominators = np.asarray(old['pose_scaling']['column_denominators'], dtype='<f8')
        if scales.shape != (384,) or denominators.shape != (7040,) or np.any(scales <= 0) or np.any(denominators <= 0):
            raise ValueError('Invalid inherited positive scales')
        if hashlib.sha256(scales.tobytes()).hexdigest() != old['pose_scaling']['group_scales_sha256']:
            raise ValueError('Changed cubic group scales')
        if hashlib.sha256(denominators.tobytes()).hexdigest() != old['pose_scaling']['column_denominators_sha256']:
            raise ValueError('Changed cubic column denominators')
        op.denominators = denominators.copy(); L = float(np.sqrt(scales.sum()))
        pilot_checkpoint = BASE/'representation/10049/real_particles-spacing-2.0.npz'
        cellop, pilot, _, check_noise = model(g, np.load(pilot_checkpoint), 24, noise=noise); pilot = cellop.expand(pilot)
        np.testing.assert_allclose(check_noise, noise); np.testing.assert_allclose(np.linalg.norm(pilot), 1., rtol=1e-12)
        moments = cell_moments_cubic(g['k'], pilot, 24, nthreads=2)
        remainder = SobolevRemainderPenalty(g['k'], g['q'], op.transfer, angle, shift, degree=3, domain='cube')
        objective = CubicDesignObjective(gram, op, remainder, moments, centers, signs, width, 2., 1., L,
            norm.isf((CONFIG['alpha_total']-CONFIG['delta'])/2), smoothing=.002)
        result.update(setup_seconds=time.perf_counter()-start, lifted_radius=L, width_fraction_field=width,
            centers_fraction_field=centers, signs=signs, field_A=g['field_A'],
            pilot_checkpoint_sha256=sha(pilot_checkpoint), remainder_gram_diagnostics=remainder.diagnostics,
            coefficient_tensor_bytes=op.polynomials.nbytes, remainder_gram_bytes=remainder.blocks.nbytes,
            preconditioner_diagnostics=gram.preconditioner_diagnostics)
        save(); progress({'stage': 'optimization', 'setup_seconds': result['setup_seconds']})
        fitted = optimize_cubic_weights(objective, w, callback=progress, checkpoint_callback=checkpoint)
        w = fitted.pop('weights'); result['fit'] = fitted
        op.set_weights(w)
        bell = [higher_order_particle_remainders(g['k'][i], g['q'][i], op.c[i], angle, shift,
            degrees=(3,), domain='cube')['records'][0]['field_remainder'] for i in range(op.n)]
        direct_remainder = float(np.sum(bell)); psd_remainder = fitted['quartic_bias']/3
        cross = cubic_residual_cross(op, w, noise, centers, signs, width); vector = cross.pop('cross_vector')
        pilot_vector = op.pair_moments(moments)
        joint = joint_density_pose_bias(fitted['density_bias']/2, fitted['pose_field_upper'], cross['norm_upper'], L,
            2., 1., min(direct_remainder, psd_remainder), float(np.linalg.norm(pilot_vector)))
        half = bias_aware_half_width_stable(fitted['noise_sd'], joint['bias_upper'], fitted['alpha_noise'])
        no_data = 2*fitted['target_norm']; pilot_target = float(cell_target_coefficients(24, centers, signs, width)@pilot)
        np.testing.assert_allclose(pilot_target, original['pilot_target'], rtol=1e-12)
        arrays = path.with_suffix('.npz')
        np.savez(arrays, weights=w, indices=g['indices'], noise_std=noise, group_scales=scales,
            column_denominators=denominators, residual_pose_cross=vector, pilot_pose_cross=pilot_vector,
            direct_quartic_remainder_per_particle=np.asarray(bell))
        result.update(bound=joint, cross=cross, direct_remainder_field=direct_remainder,
            psd_remainder_field=psd_remainder, refined_half_width=half, no_data_half_width=no_data,
            selected_relative_half_width=min(1., half/no_data), uses_no_data=bool(half >= no_data),
            pilot_target=pilot_target, arrays_sha256=sha(arrays), reference_map_sha256=sha(ROOT/'data/uncertainty/references/emd_6487.map'),
            interval_selection_scope='Original quadratic, fixed cubic and optimized cubic are separate alternatives, not an unadjusted selected joint interval.')
        save(); progress({'stage': 'nonlinear_reference_checks', 'relative_width': result['selected_relative_half_width'],
                          'surrogate_gap': fitted['relative_sum_gap']})
        nominal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
        rho = VoxelReference.from_mrc(ROOT/'data/uncertainty/references/emd_6487.map', box=64).volume.ravel()
        rho /= np.linalg.norm(rho); truth = float(cell_target_coefficients(64, centers, signs, width)@rho)
        rng = np.random.default_rng(610281+10049)
        for scenario in ['nominal', 'coherent_x', 'random_boundary']:
            poses = np.zeros((op.n, 5))
            if scenario == 'coherent_x': poses[:, 0] = 1.
            if scenario == 'random_boundary':
                poses = rng.normal(size=poses.shape); poses /= np.linalg.norm(poses, axis=1)[:, None]
            signal = pose_cell_forward(g['k'], g['q'], g['ctf'], rho, 64, noise, poses, angle, shift)
            expected = float(pilot_target+w@(signal-nominal))
            check = reference_interval_summary(truth, expected, fitted['noise_sd'], pilot_target, half, no_data)
            result['reference_checks'].append(dict(scenario=scenario, raw_expected_center=expected, **check)); save()
        result.update(complete=True, seconds=time.perf_counter()-start,
            minimum_reference_power=min(r['correct_sign_probability'] for r in result['reference_checks']),
            reference_coverage_failure=any(r['analytic_coverage'] < 1-fitted['alpha_noise']-1e-8 for r in result['reference_checks']),
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024))
        save(); print('COMPLETE', result['selected_relative_half_width'], result['minimum_reference_power'], flush=True)
    except Exception as exc:
        result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
