#!/usr/bin/env python3
"""Fresh-design validation of the continuous audit under a locked protocol.

A smoke option exercises the same code on old development geometry only. The
ordinary mode requires a committed protocol and source/data lock manifest.
"""
import argparse
import csv
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm, beta
from fourier_splats.physics import ctf
from fourier_splats.uq_data import particle_geometry, half_plane_frequencies, VoxelReference
from fourier_splats.uq_continuous import (continuous_certificate, continuous_residual_norm, cell_forward,
                                        cell_adjoint, cell_target_coefficients)
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_continuous_pose import continuous_pose_audit, pose_cell_forward, perturbed_geometry
from fourier_splats.uncertainty import bias_aware_half_width
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model
ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT/'research/uncertainty/confirmation/continuous-v1'
MAPS = {'10028': '2660', '10049': '6487', '10076': '8434'}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as file:
        for block in iter(lambda: file.read(16*1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def procedural_density(seed, box, kind):
    """Independent cell-defined shapes; no Fourier-Gaussian fitting dictionary."""
    rng = np.random.default_rng(seed); axis = (np.arange(box)+.5)/box-.5
    z, y, x = np.meshgrid(axis, axis, axis, indexing='ij'); xyz = np.stack([x, y, z], axis=-1)
    density = np.zeros((box,)*3)
    if kind == 'ellipsoids':
        for _ in range(12):
            center = rng.uniform(-.24, .24, 3); radii = rng.uniform(.04, .16, 3)
            density += rng.uniform(.2, 1.)*(np.sum(((xyz-center)/radii)**2, axis=-1) <= 1)
    elif kind == 'modulated_shell':
        center = rng.uniform(-.07, .07, 3); scale = rng.uniform(.11, .23, 3)
        distance = np.sum(((xyz-center)/scale)**2, axis=-1)
        density = np.exp(-((distance-1)/.5)**2)
        for _ in range(5):
            frequency = rng.integers(3, 13, 3); phase = rng.uniform(0, 2*np.pi)
            density *= 1+.12*np.cos(2*np.pi*np.sum(xyz*frequency, axis=-1)+phase)
    elif kind == 'signed_multiscale':
        envelope = np.exp(-np.sum(xyz**2, axis=-1)/.12)
        for _ in range(12):
            frequency = rng.integers(1, 15, 3); phase = rng.uniform(0, 2*np.pi)
            density += rng.normal()*np.cos(2*np.pi*np.sum(xyz*frequency, axis=-1)+phase)
        density *= envelope
    else:
        raise ValueError('Unknown generator')
    density = density.ravel(); density /= np.linalg.norm(density)
    return density


def fresh_geometry(dataset, replicate, config):
    data = ROOT/'data/uncertainty/confirmation/prediction-v1'/dataset
    metadata = np.load(data/'metadata.npz'); manifest = json.loads((data/'manifest.json').read_text())
    count = config['particles']; indices = np.random.default_rng(config['geometry_seed']+int(dataset)).permutation(4096)
    indices = indices[replicate*count:(replicate+1)*count]
    q = half_plane_frequencies(config['frequency_radius']); plane = np.pad(q, ((0, 0), (0, 1)))
    k = np.einsum('qi,nij->nqj', plane, metadata['rotations'][indices])
    field = manifest['raw_box']*manifest['raw_pixel_size_A']
    return {'indices': indices, 'q': np.broadcast_to(q, (count, len(q), 2)).copy(), 'k': k,
            'ctf': ctf(q/field, metadata['ctf'][indices]).astype(float), 'field_A': field}


def coverage_record(half, actual_bias, sd, draws):
    if sd > 0:
        analytic = float(norm.cdf((half-actual_bias)/sd)-norm.cdf((-half-actual_bias)/sd))
        count = int(np.sum(abs(actual_bias+sd*draws) <= half))
    else:
        analytic = float(abs(actual_bias) <= half+1e-12); count = int(analytic*len(draws))
    n = len(draws)
    interval = [0. if count == 0 else float(beta.ppf(.025, count, n-count+1)),
                1. if count == n else float(beta.ppf(.975, count+1, n-count))]
    return {'half_width': float(half), 'actual_bias': float(actual_bias), 'actual_noise_sd': float(sd),
            'analytic_coverage': analytic, 'covered': count, 'replications': n, 'monte_carlo_coverage': count/n,
            'binomial_95': interval}


def run(config, dataset, replicate, smoke, snapshot, lock_hash):
    start = time.perf_counter()
    out = ROOT/('results/uncertainty/development/continuous-confirmation-smoke' if smoke else 'results/uncertainty/confirmation/continuous-v1')
    out.mkdir(parents=True, exist_ok=True); path = out/f'{dataset}-replicate{replicate}.json'
    if path.exists():
        old = json.loads(path.read_text())
        if old['protocol_lock_sha256'] != lock_hash:
            raise AssertionError('Cannot resume under a different protocol')
        if old.get('complete'):
            print(dataset, replicate, 'already complete', flush=True); return
    g = particle_geometry(ROOT, dataset, 'inference_half0', radius=5, count=config['particles'], seed=900001) if smoke else fresh_geometry(dataset, replicate, config)
    checkpoint_path = ROOT/'results/uncertainty/development/representation'/dataset/'real_particles-spacing-2.0.npz'
    checkpoint = np.load(checkpoint_path); op, pilot, _, noise = model(g, checkpoint, 24)
    pilot = op.expand(pilot); nominal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
    B = config['density_radius']; fine = config['generator_box']; P = float(np.linalg.norm(pilot))
    gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=config['quadrature_order'], preconditioner_rank=config['preconditioner_rank'])
    reference = VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map', box=fine).volume.ravel()
    reference /= np.linalg.norm(reference)
    generators = [('deposited_map', reference)]
    for offset, kind in enumerate(['ellipsoids', 'modulated_shell', 'signed_multiscale']):
        seed = config['signal_seed']+int(dataset)*100+replicate*10+offset
        generators.append((kind, procedural_density(seed, fine, kind)))
    result = {'stage': 'development smoke only' if smoke else 'frozen continuous-v1 simulated coverage on additional-exposure acquisition geometries',
              'dataset': dataset, 'replicate': replicate, 'config': config, 'source_snapshot': snapshot,
              'protocol_lock_sha256': lock_hash, 'complete': False, 'geometry_indices': g['indices'].tolist(),
              'geometry_sha256': hashlib.sha256(g['k'].tobytes()+g['ctf'].tobytes()).hexdigest(),
              'noise_std_physical': float(noise), 'generator_hashes': {name: hashlib.sha256(rho.tobytes()).hexdigest() for name, rho in generators},
              'records': []}
    rng = np.random.default_rng(config['noise_seed']+int(dataset)*100+replicate)
    for width in config['widths']:
        for name in config['targets']:
            centers = [[0, 0, 0]] if name == 'center' else [[0, 0, .08], [0, 0, -.08]]; signs = [1] if name == 'center' else [1, -1]
            fit = continuous_certificate(gram, centers, signs, width, B, maxiter=100, rtol=.005)
            w = fit.pop('weights'); sd = float(np.linalg.norm(w)); no_data = B*fit['target_norm']
            exact_nominal = continuous_residual_norm(g['k'], g['ctf'], w, noise, centers, signs, width)
            if abs(exact_nominal['residual_norm']-fit['bias']/B) > 1e-5:
                raise AssertionError('Continuous primal disagreement')
            target_pilot = float(cell_target_coefficients(24, centers, signs, width)@pilot)
            fine_target = cell_target_coefficients(fine, centers, signs, width)
            cell_h = cell_target_coefficients(24, centers, signs, width)-cell_adjoint(g['k'], g['ctf'], w, 24, noise)
            baseline_widths = {'noise_only': float(norm.ppf(.975)*sd), 'cell24_fixed_pose': float(bias_aware_half_width(sd, B*np.linalg.norm(cell_h))),
                               'continuous_fixed_pose': float(fit['half_width'])}
            a, target_norm2 = gram.target(centers, signs, width)
            nominal_hnorm = np.sqrt(max(0., exact_nominal['residual_norm2_unpadded']))
            nominal_pairing = target_norm2-float(w@a)
            for degrees in config['angles']:
                angle = np.deg2rad(degrees); shift = config['shift_field']
                audit = continuous_pose_audit(g['k'], g['q'], g['ctf'], w, noise, angle, shift, B, P, fit['bias']/B,
                                             order=config['pose_quadrature_order'])
                full_half = float(bias_aware_half_width(sd, audit['total_bias']))
                widths = dict(baseline_widths, continuous_pose=full_half, no_data=float(no_data))
                row = {'target': name, 'width': width, 'angle_degrees': degrees, 'fit': fit, 'audit': audit,
                       'half_widths_before_fallback': widths, 'no_data_half_width': float(no_data), 'cases': []}
                # Pose draws use a stream separate from scalar measurement noise.
                pose_rng = np.random.default_rng(config['pose_seed']+int(dataset)*100+replicate*10+int(degrees*10))
                random_pose = pose_rng.normal(size=(len(g['k']), 5)); random_pose /= np.linalg.norm(random_pose, axis=1)[:, None]
                coherent = np.zeros_like(random_pose); coherent[:, 0] = 1
                scenarios = [('nominal', np.zeros_like(random_pose), 1., True), ('random_boundary', random_pose, 1., True),
                             ('coherent_boundary', coherent, 1., True), ('pose_exceedance4', 4*coherent, 1., False),
                             ('noise_scale2', random_pose, 2., False)]
                signal_cases = []
                for label, pose, scale, in_class in scenarios:
                    for generator, rho in generators:
                        observed = pose_cell_forward(g['k'], g['q'], g['ctf'], rho, fine, noise, pose, angle, shift)
                        truth = float(fine_target@rho); expected = float(target_pilot+w@(observed-nominal))
                        signal_cases.append((generator, label, truth, expected-truth, scale, in_class))
                # Full-L2 nominal boundary signs, not a cell-projected adversary.
                for sign in [-1, 1]:
                    truth = target_pilot+sign*B*nominal_pairing/nominal_hnorm
                    signal_cases.append((f'continuous_nominal_boundary_{sign}', 'nominal', truth, -sign*B*nominal_hnorm, 1., True))
                rotated, translations = perturbed_geometry(g['k'], g['q'], coherent, angle, shift)
                moving = continuous_residual_norm(rotated, g['ctf']*np.exp(2j*np.pi*translations), w, noise, centers, signs, width)
                moving_norm = np.sqrt(max(0., moving['residual_norm2_unpadded']))
                pilot_signal = pose_cell_forward(g['k'], g['q'], g['ctf'], pilot, 24, noise, coherent, angle, shift)
                pilot_bias = float(w@(pilot_signal-nominal))
                for sign in [-1, 1]:
                    truth = target_pilot+sign*B*(moving['cross_term']-target_norm2)/moving_norm
                    signal_cases.append((f'continuous_pose_boundary_{sign}', 'coherent_boundary', truth, pilot_bias+sign*B*moving_norm, 1., True))
                for generator, scenario, truth, actual_bias, noise_scale, in_class in signal_cases:
                    draws = rng.standard_normal(config['noise_replicates']); methods = {}
                    for method, half in widths.items():
                        fallback = method == 'no_data' or half >= no_data
                        bias = target_pilot-truth if fallback else actual_bias
                        current_sd = 0. if fallback else sd*noise_scale
                        methods[method] = coverage_record(min(half, no_data), bias, current_sd, draws)
                        methods[method]['uses_no_data'] = fallback
                        methods[method]['correct_sign_probability'] = (float(norm.cdf((np.sign(truth)*(truth+bias)-min(half, no_data))/current_sd))
                            if current_sd > 0 else float(np.sign(truth)*(truth+bias) > min(half, no_data)))
                    if in_class and methods['continuous_pose']['analytic_coverage'] < .95-1e-8:
                        raise AssertionError('In-class continuous pose undercoverage')
                    row['cases'].append({'generator': generator, 'scenario': scenario, 'within_declared_model': in_class,
                                         'true_target': float(truth), 'noise_scale': noise_scale, 'methods': methods})
                result['records'].append(row); result['seconds'] = time.perf_counter()-start
                path.write_text(json.dumps(result, indent=2)+'\n')
                print(dataset, replicate, name, width, degrees, 'cases', len(row['cases']), 'relative width', min(full_half/no_data, 1.), flush=True)
            np.savez(out/f'{dataset}-replicate{replicate}-{name}-{width}-weights.npz', weights=w)
    result['complete'] = True; result['seconds'] = time.perf_counter()-start
    path.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--smoke', action='store_true'); args = parser.parse_args()
    snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py'])
    if args.smoke:
        config = {'particles': 8, 'frequency_radius': 5., 'density_radius': 2., 'generator_box': 24,
                  'quadrature_order': 32, 'pose_quadrature_order': 24, 'preconditioner_rank': 128,
                  'widths': [.07], 'targets': ['center'], 'angles': [.1], 'shift_field': .01/24,
                  'signal_seed': 900011, 'pose_seed': 900021, 'noise_seed': 900031, 'noise_replicates': 100}
        run(config, '10028', 0, True, snapshot, 'development-smoke-'+sha(Path(__file__)))
    else:
        lock_path = PROTOCOL/'locked-study.json'; locked = json.loads(lock_path.read_text())
        for path, expected in locked['file_hashes'].items():
            if sha(ROOT/path) != expected:
                raise AssertionError(f'Frozen study dependency drift: {path}')
        for dataset in locked['datasets']:
            for replicate in range(locked['replicates']):
                run(locked['config'], dataset, replicate, False, snapshot, sha(lock_path))
