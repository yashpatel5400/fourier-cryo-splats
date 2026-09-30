#!/usr/bin/env python3
"""Fit every locked pilot-only target under the specified continuous class."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import finufft
import numpy as np
from fourier_splats.uq_continuous import continuous_certificate, cell_forward, cell_target_coefficients
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_data import particle_geometry, VoxelReference
from fourier_splats.uq_intervals import bias_aware_half_width_stable, reference_interval_summary
from fourier_splats.uq_pilot_targets import read_target_lock, verify_pilot_scores, LOCK_SHA256
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
LOCK = ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json'
MAPS = {'10028': '2660', '10049': '6487', '10076': '8434'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(args, locked, snapshot):
    dataset = locked['dataset']; out = BASE/args.output; out.mkdir(parents=True, exist_ok=True)
    requested = args.features.split(',') if args.features else [f['name'] for f in locked['features']]
    if set(requested)-{f['name'] for f in locked['features']}:
        raise ValueError('Target not in lock')
    features = []
    for feature in locked['features']:
        if feature['name'] not in requested:
            continue
        path = out/f"{dataset}-{feature['name']}.json"
        if path.exists():
            previous = json.loads(path.read_text())
            if not args.resume or not previous.get('complete') or previous.get('error'):
                raise RuntimeError(f'Preserve existing outcome: {path}')
            if previous['target_lock_sha256'] != LOCK_SHA256:
                raise ValueError('Existing outcome uses a different target lock')
            print('REUSE', path.name, flush=True)
        else:
            features.append(feature)
    if not features:
        return
    setup = time.perf_counter()
    g = particle_geometry(ROOT, dataset, 'inference_half0', radius=12, count=128, seed=609315)
    np.testing.assert_allclose(g['field_A'], locked['field_A'], rtol=1e-14)
    checkpoint = BASE/'representation'/dataset/'real_particles-spacing-2.0.npz'
    if sha(checkpoint) != locked['pilot_checkpoint_sha256']:
        raise ValueError('Pilot checkpoint changed after target selection')
    noise_path = BASE/'continuous-quadrature-optimized'/f'{dataset}-center-0.07-weights.npz'
    noise = float(np.load(noise_path)['noise_std'])
    operator, pilot, _, _ = model(g, np.load(checkpoint), 24, noise=noise)
    pilot = operator.expand(pilot); verify_pilot_scores(locked, pilot)
    gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=80, preconditioner_rank=1024)
    gram.nthreads = args.threads
    nominal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
    setup_seconds = time.perf_counter()-setup
    alpha_total = .05/12; delta = 1e-6/12; alpha_noise = alpha_total-delta
    for feature in features:
        name = feature['name']; path = out/f'{dataset}-{name}.json'; start = time.perf_counter()
        centers = [feature['center_fraction_field']]; signs = [1]; width = locked['width_fraction_field']
        result = {'complete': False, 'dataset': dataset,
                  'stage': 'Exploratory pilot-selected targets; fixed poses/CTFs and prescribed white Gaussian noise',
                  'class': 'radius 2 around a unit-L2 constant-cell pilot; independent unit-L2 64-cell reference generator',
                  'density_radius': 2., 'pilot_norm_bound': 1., 'target_lock_sha256': LOCK_SHA256,
                  'config': dict(vars(args), particles=128, frequency_radius=12., seed=609315, quadrature_order=80,
                                 preconditioner_rank=1024, alpha_total=alpha_total, alpha_noise=alpha_noise, delta=delta),
                  'source_snapshot': snapshot, 'pilot_checkpoint_sha256': sha(checkpoint),
                  'noise_source_sha256': sha(noise_path), 'noise_scope': 'Prescribed original radius-five pilot-SNR scale, not experimental calibration',
                  'noise_std': noise, 'field_A': g['field_A'], 'band_endpoint_A': g['field_A']/12,
                  'setup_seconds_shared_per_dataset': setup_seconds, 'targets': []}
        def save():
            path.write_text(json.dumps(result, indent=2)+'\n')
        save()
        try:
            fit = continuous_certificate(gram, centers, signs, width, 2., alpha=alpha_noise, maxiter=100,
                rtol=.005, callback=lambda r: print(dataset, name, r, flush=True))
            w = fit.pop('weights')
            fit['half_width'] = bias_aware_half_width_stable(fit['noise_sd'], fit['bias'], alpha_noise)
            wp = out/f'{dataset}-{name}-{width}-weights.npz'
            np.savez(wp, weights=w, indices=g['indices'], noise_std=noise)
            row = {'target': name, 'centers_fraction_field': centers, 'signs': signs, 'width_fraction_field': width,
                   'width_A': 20., 'pilot_expected_feature': feature['pilot_expected_feature'], 'fit': fit,
                   'source_weights_sha256': sha(wp), 'solve_seconds': time.perf_counter()-start,
                   'selected_relative_half_width': min(1., fit['half_width']/(2*fit['target_norm']))}
            result['targets'].append(row); save()
            # First reference read for this feature occurs only after its weights are saved.
            rp = ROOT/f'data/uncertainty/references/emd_{MAPS[dataset]}.map'
            reference = VoxelReference.from_mrc(rp, box=64).volume.ravel(); reference /= np.linalg.norm(reference)
            truth = float(cell_target_coefficients(64, centers, signs, width)@reference)
            expected = float(feature['pilot_expected_feature']+w@(cell_forward(g['k'], g['ctf'], reference, 64, noise)-nominal))
            check = reference_interval_summary(truth, expected, fit['noise_sd'], feature['pilot_expected_feature'],
                                               fit['half_width'], 2*fit['target_norm'])
            row.update(reference_check=check, raw_expected_center=expected, reference_sha256=sha(rp))
            result.update(complete=True, seconds=time.perf_counter()-start,
                          reference_coverage_failure=bool(check['analytic_coverage'] < 1-alpha_noise-1e-8))
            save(); print('DONE', dataset, name, row['selected_relative_half_width'], check['correct_sign_probability'], flush=True)
        except Exception as exc:
            result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--datasets', default='10028,10049,10076'); parser.add_argument('--features')
    parser.add_argument('--threads', type=int, default=2); parser.add_argument('--resume', action='store_true')
    parser.add_argument('--output', default='pilot-selected-fixed-v1'); args = parser.parse_args()
    if args.threads < 1 or (args.threads > 1 and ('tmp/nufft-openmp' not in finufft.__file__ or 'torch' in sys.modules)):
        raise RuntimeError('Use isolated CPU-only FINUFFT for a threaded Mac run')
    lock = read_target_lock(LOCK)
    if set(args.datasets.split(','))-{r['dataset'] for r in lock['datasets']}:
        raise ValueError('Dataset not in target lock')
    snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py',
        'research/uncertainty/pilot-selected-targets-v1/EXPERIMENT.md', str(LOCK.relative_to(ROOT))])
    for row in lock['datasets']:
        if row['dataset'] in args.datasets.split(','):
            run(args, row, snapshot)


if __name__ == '__main__':
    main()
