#!/usr/bin/env python3
"""Declared joint density/pose cubic design; retains every predecessor."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time
import finufft
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_cubic_optimization import CubicDesignObjective
from fourier_splats.uq_cubic_subspace import orthonormal_basis
from fourier_splats.uq_joint_cubic_design import design_joint_cubic, audit_joint_cubic
from fourier_splats.uq_cubic_design import DifferentiableCubicPoseFieldOperator
from fourier_splats.uq_cubic_pose import cell_moments_cubic
from fourier_splats.uq_sobolev_penalty import SobolevRemainderPenalty
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_continuous import cell_forward, cell_target_coefficients
from fourier_splats.uq_continuous_pose import pose_cell_forward
from fourier_splats.uq_data import particle_geometry, VoxelReference
from fourier_splats.uq_intervals import reference_interval_summary
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
SOURCE = BASE/'cubic-pose-probe/10049-pilot_region_1-1.json'
PROTOCOL = ROOT/'research/uncertainty/JOINT-CUBIC-DESIGN-PROTOCOL.md'
CONFIG = {'dataset': '10049', 'target': 'pilot_region_1', 'particles': 128, 'radius': 12,
          'width_A': 20., 'angle_degrees': 1., 'shift_A': .5, 'density_radius': 2., 'pilot_norm': 1.,
          'density_order': 80, 'pose_order': 64, 'preconditioner_rank': 1024, 'threads': 2,
          'evaluations': 12, 'krylov_steps': 24, 'audit_krylov_steps': 40, 'design_seconds': 5400, 'total_alarm_seconds': 10800, 'separation_tolerance': .001,
          'optimization_seed': 953011, 'certificate_seed': 953001, 'probes': 4, 'power_iterations': 40,
          'alpha_total': .05/12, 'delta': 1e-6/12, 'basis_relative_svd_tolerance': 1e-12}



def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seed_inventory():
    matches = []; count = 0
    def scan(value, keys, name):
        if isinstance(value, dict):
            for key, child in value.items():
                if 'seed' in key and isinstance(child, int) and child in [CONFIG['certificate_seed'], CONFIG['optimization_seed']]:
                    matches.append({'path': name, 'keys': keys+[key]})
                scan(child, keys+[key], name)
        elif isinstance(value, list):
            for j, child in enumerate(value):
                scan(child, keys+[j], name)
    for path in sorted(BASE.rglob('*.json')):
        count += 1; scan(json.loads(path.read_text()), [], str(path.relative_to(ROOT)))
    if matches:
        raise ValueError(f'Certificate seed already recorded: {matches}')
    return {'json_files_examined': count, 'matches': matches, 'seeds': [CONFIG['certificate_seed'], CONFIG['optimization_seed']],
            'scope': 'All existing development JSON seed-valued fields before creating this run; not text substring matches.'}


def main():
    if 'tmp/nufft-openmp' not in finufft.__file__ or 'torch' in sys.modules:
        raise RuntimeError('Isolated threaded CPU-only FINUFFT required')
    out = BASE/'joint-cubic-design-probe'; path = out/'10049-pilot_region_1-1.json'
    if path.exists():
        raise RuntimeError('Preserve existing attempt, including failures')
    parser = argparse.ArgumentParser(); parser.add_argument('--wait-hours', type=float, default=0.)
    args = parser.parse_args(); deadline = time.monotonic()+3600*args.wait_hours
    predecessor_path = BASE/'cubic-weight-probe/10049-pilot_region_1-1.json'
    predecessor = json.loads(predecessor_path.read_text())
    if not predecessor.get('complete') or predecessor.get('error'):
        raise RuntimeError('Completed original cubic fit required')
    print('WAITING_FOR_ENRICHMENT_PROCESS_TERMINATION', flush=True)
    coordinate_path = BASE/'cubic-enrichment-probe/10049-pilot_region_1-1.json'
    while True:
        try:
            scheduling = json.loads(coordinate_path.read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            scheduling = {}
        if scheduling.get('complete') or scheduling.get('error'):
            break
        if time.monotonic() >= deadline:
            raise TimeoutError('Enrichment run has not terminated; do not oversubscribe heavy fits')
        time.sleep(15)
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
    extra = ['scripts/audit_uq_grid_refinement.py', 'research/uncertainty/JOINT-CUBIC-DESIGN-PROTOCOL.md',
             'tests/test_uq_joint_cubic_design.py', 'tests/test_uq_trust_region.py',
             'research/uncertainty/JOINT-TRUST-REGION-THEORY.md', 'research/uncertainty/CUBIC-WEIGHT-RESULTS.md']
    # No empirical run may silently execute unpublished changed implementation.
    for p in [Path(__file__), *sorted((ROOT/'src/fourier_splats').glob('*.py')), *[ROOT/f for f in extra]]:
        relative = str(p.relative_to(ROOT))
        if subprocess.check_output(['git', 'show', f'HEAD:{relative}'], cwd=ROOT) != p.read_bytes():
            raise ValueError(f'Commit declared source before empirical fitting: {relative}')
    out.mkdir(parents=True, exist_ok=True)
    result = {'complete': False, 'config': CONFIG, 'dataset': '10049', 'target': 'pilot_region_1',
        'scope': 'Joint density/pose guide from the unchanged thirteen-column span and original cubic fit; same class/noise/pose assumptions, not fresh experimental calibration.',
        'predecessor_case': str(predecessor_path.relative_to(ROOT)), 'predecessor_sha256': sha(predecessor_path),
        'predecessor_outcome': {k: predecessor[k] for k in ['selected_relative_half_width', 'minimum_reference_power']},
        'predecessor_relative_gap': predecessor['fit']['relative_sum_gap'],
        'source_case': str(SOURCE.relative_to(ROOT)), 'source_case_sha256': sha(SOURCE),
        'source_fit': str(source_fit.relative_to(ROOT)), 'source_fit_sha256': sha(source_fit),
        'source_weights_sha256': sha(wp), 'protocol_sha256': sha(PROTOCOL), 'seed_inventory': inventory,
        'source_snapshot': source_snapshot(ROOT, Path(__file__), extra),
        'checkpoints': [], 'reference_checks': [], 'pose_scaling': old['pose_scaling']}
    start = time.perf_counter()
    def timeout_handler(signum, frame):
        raise TimeoutError('Declared three-hour total alarm reached')
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(CONFIG['total_alarm_seconds'])
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
            norm.isf((CONFIG['alpha_total']-CONFIG['delta'])/2), smoothing=0.)
        result.update(setup_seconds=time.perf_counter()-start, lifted_radius=L, width_fraction_field=width,
            centers_fraction_field=centers, signs=signs, field_A=g['field_A'],
            pilot_checkpoint_sha256=sha(pilot_checkpoint), remainder_gram_diagnostics=remainder.diagnostics,
            coefficient_tensor_bytes=op.polynomials.nbytes, remainder_gram_bytes=remainder.blocks.nbytes,
            preconditioner_diagnostics=gram.preconditioner_diagnostics)
        save(); progress({'stage': 'optimization', 'setup_seconds': result['setup_seconds']})
        stamp = time.perf_counter()
        original_array = predecessor_path.with_suffix('.npz')
        if sha(original_array) != predecessor['arrays_sha256']:
            raise ValueError('Original optimized weights changed')
        initial_weights = np.load(original_array)['weights']
        shell = np.ceil(np.linalg.norm(g['q'], axis=-1)-1e-12).astype(int)
        if shell.min() < 1 or shell.max() > 12:
            raise ValueError('Unexpected Fourier shell labels')
        packed = w.reshape(op.n, 2, op.nq)
        columns = [np.where((shell == j)[:,None,:], packed, 0.).ravel() for j in range(1,13)]
        columns.append(initial_weights)
        basis, basis_record = orthonormal_basis(np.column_stack(columns), CONFIG['basis_relative_svd_tolerance'])
        result.update(basis_seconds=time.perf_counter()-stamp, basis_diagnostics=basis_record,
            original_optimized_array_sha256=sha(original_array))
        save(); progress({'stage': 'reduced_basis_ready', 'seconds': result['basis_seconds'], 'rank': basis.shape[1]})
        designed = design_joint_cubic(objective, basis, initial_weights, centers, signs, width,
            optimization_seed=CONFIG['optimization_seed'], evaluations=CONFIG['evaluations'],
            krylov_steps=CONFIG['krylov_steps'], design_seconds=CONFIG['design_seconds'],
            separation_tolerance=CONFIG['separation_tolerance'], callback=progress, checkpoint_callback=checkpoint)
        w = designed.pop('weights'); result['design'] = designed
        save(); progress({'stage': 'fresh_joint_audit', 'selected_evaluation': designed['selected_evaluation']})
        fitted = audit_joint_cubic(objective, w, centers, signs, width,
            certificate_seed=CONFIG['certificate_seed'], delta=CONFIG['delta'], alpha=CONFIG['alpha_total'],
            probes=CONFIG['probes'], power_iterations=CONFIG['power_iterations'],
            krylov_steps=CONFIG['audit_krylov_steps'], callback=progress)
        vector = fitted.pop('residual_cross_vector'); result['fit'] = fitted
        half = fitted['half_width']
        no_data = 2*fitted['target_norm']; pilot_target = float(cell_target_coefficients(24, centers, signs, width)@pilot)
        np.testing.assert_allclose(pilot_target, original['pilot_target'], rtol=1e-12)
        arrays = path.with_suffix('.npz')
        np.savez(arrays, weights=w, reduced_basis=basis, indices=g['indices'], noise_std=noise, group_scales=scales,
            column_denominators=denominators, residual_pose_cross=vector)
        result.update(refined_half_width=half, no_data_half_width=no_data,
            selected_relative_half_width=min(1., half/no_data), uses_no_data=bool(half >= no_data),
            pilot_target=pilot_target, arrays_sha256=sha(arrays), reference_map_sha256=sha(ROOT/'data/uncertainty/references/emd_6487.map'),
            interval_selection_scope='Joint and all preceding cubic runs remain separate alternatives; no unadjusted minimum across procedures is claimed.')
        save(); progress({'stage': 'nonlinear_reference_checks', 'relative_width': result['selected_relative_half_width'],
            'restricted_guide_gap': designed['restricted_guide_gap']})
        nominal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
        original_rho = VoxelReference.from_mrc(ROOT/'data/uncertainty/references/emd_6487.map', box=64).volume.ravel()
        registration_path = BASE/'reference-registration-v1/10049.json'
        registration = json.loads(registration_path.read_text())
        registration_arrays = ROOT/registration['arrays_file']
        if sha(registration_arrays) != registration['arrays_sha256']:
            raise ValueError('Registered reference changed')
        registered_rho = np.load(registration_arrays)['reference_registered'].ravel().copy()
        result['registered_reference_source_sha256'] = sha(registration_path)
        result['registered_reference_arrays_sha256'] = sha(registration_arrays)
        for frame, rho in [('original', original_rho), ('registered', registered_rho)]:
            rho /= np.linalg.norm(rho); truth = float(cell_target_coefficients(64, centers, signs, width)@rho)
            rng = np.random.default_rng(610281+10049)
            for scenario in ['nominal', 'coherent_x', 'random_boundary']:
                poses = np.zeros((op.n, 5))
                if scenario == 'coherent_x': poses[:, 0] = 1.
                if scenario == 'random_boundary':
                    poses = rng.normal(size=poses.shape); poses /= np.linalg.norm(poses, axis=1)[:, None]
                signal_values = pose_cell_forward(g['k'], g['q'], g['ctf'], rho, 64, noise, poses, angle, shift)
                expected = float(pilot_target+w@(signal_values-nominal))
                check = reference_interval_summary(truth, expected, fitted['noise_sd'], pilot_target, half, no_data)
                result['reference_checks'].append(dict(frame=frame, scenario=scenario, raw_expected_center=expected, **check)); save()
        result['minimum_reference_power_by_frame'] = {frame: min(r['correct_sign_probability'] for r in result['reference_checks'] if r['frame']==frame)
            for frame in ['original','registered']}
        result.update(complete=True, seconds=time.perf_counter()-start,
            minimum_reference_power=min(r['correct_sign_probability'] for r in result['reference_checks']),
            reference_coverage_failure=any(r['analytic_coverage'] < 1-fitted['alpha_noise']-1e-8 for r in result['reference_checks']),
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024))
        signal.alarm(0)
        save(); print('COMPLETE', result['selected_relative_half_width'], result['minimum_reference_power'], flush=True)
    except Exception as exc:
        signal.alarm(0)
        result.update(complete=True, stage='failed', error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
