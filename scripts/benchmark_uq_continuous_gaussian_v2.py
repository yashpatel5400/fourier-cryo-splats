#!/usr/bin/env python3
"""Review-2 matched continuous-prior comparison; all declared outcomes retained."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from scipy.stats import norm
from fourier_splats.uq_continuous import cell_forward, cell_target_coefficients
from fourier_splats.uq_continuous_gaussian import continuous_gaussian, cell_pose_jacobian
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_cached_quadrature import CachedQuadratureGram
from fourier_splats.uq_continuous_pose import pose_cell_forward
from fourier_splats.uq_fourier_pose_baseline import GaussianNuisanceWhitening
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_pilot_targets import read_target_lock, verify_pilot_scores, LOCK_SHA256
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
PROTOCOL = ROOT/'research/uncertainty/CONTINUOUS-GAUSSIAN-V2-NUMERICAL-AMENDMENT.md'
LOCK = ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(locked, args, snapshot):
    ds = locked['dataset']; out = BASE/args.output; out.mkdir(parents=True, exist_ok=True)
    path = out/f'{ds}.json'
    if path.exists():
        raise RuntimeError(f'Preserve existing outcome: {path}')
    start = time.perf_counter()
    result = dict(complete=False, dataset=ds, stage='Post-review-2 declared development comparison',
        scope='Matched continuous operator and directional prior scale. Conditional prescribed-noise simulations; no end-to-end or experimental coverage.',
        config=dict(particles=128, frequency_radius=12, geometry_seed=609315, prior_directional_sds=[2., 1.],
            density_radius=2., alpha=.05/12, quadrature_order=80, preconditioner_rank=8192,
            pose_coordinate_sd_degrees=1., shift_coordinate_sd_A=.5, threads=args.threads,
            cg_rtol=1e-10, maxiter=2000),
        source_snapshot=snapshot, target_lock_sha256=LOCK_SHA256, input_hashes={}, records=[])
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    def read(p):
        result['input_hashes'][str(p.relative_to(ROOT))] = sha(p)
        return json.loads(p.read_text())
    def array(p, expected=None):
        digest = sha(p)
        if expected and digest != expected:
            raise ValueError(f'Changed source array: {p}')
        result['input_hashes'][str(p.relative_to(ROOT))] = digest
        return np.load(p)
    save()
    try:
        for input_path in [ROOT/f'data/{ds}/metadata.npz', ROOT/f'data/{ds}/manifest.json', ROOT/f'research/uncertainty/splits/{ds}.csv']:
            result['input_hashes'][str(input_path.relative_to(ROOT))] = sha(input_path)
        g = particle_geometry(ROOT, ds, 'inference_half0', radius=12, count=128, seed=609315)
        noise_path = BASE/'continuous-quadrature-optimized'/f'{ds}-center-0.07-weights.npz'
        noise = float(array(noise_path)['noise_std'])
        cp = BASE/f'representation/{ds}/real_particles-spacing-2.0.npz'
        op, pilot, _, _ = model(g, array(cp, locked['pilot_checkpoint_sha256']), 24, noise=noise)
        pilot = op.expand(pilot); verify_pilot_scores(locked, pilot)
        reg = read(BASE/f'reference-registration-v1/{ds}.json')
        maps = array(ROOT/reg['arrays_file'], reg['arrays_sha256'])
        references = {}
        for frame in ['original', 'registered']:
            rho = maps['reference_'+frame].ravel().copy()
            references[frame] = rho/np.linalg.norm(rho)
        gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=80, preconditioner_rank=8192)
        gram = CachedQuadratureGram(gram, nthreads=args.threads)
        nominal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
        jacobian = cell_pose_jacobian(g['k'], g['q'], g['ctf'], pilot, 24, noise)
        scaled = jacobian*np.array([np.deg2rad(1.)]*3+[.5/g['field_A']]*2)
        nuisance = GaussianNuisanceWhitening(scaled)
        ap = out/f'{ds}-design.npz'
        np.savez_compressed(ap, scaled_jacobian=scaled, indices=g['indices'], noise_std=noise)
        signals = []
        # Same prescribed directions as the previous trilinear baseline.
        for frame, rho in references.items():
            signals.append((frame, 'nominal', cell_forward(g['k'], g['ctf'], rho, 64, noise)))
            for angle in [0., 1., 2.]:
                rng = np.random.default_rng(610281+int(ds))
                for direction in ['coherent_x', 'random_boundary']:
                    u = np.zeros((128, 5))
                    if direction == 'coherent_x':
                        u[:, 0] = 1.
                    else:
                        u = rng.normal(size=u.shape); u /= np.linalg.norm(u, axis=1)[:, None]
                    signal = pose_cell_forward(g['k'], g['q'], g['ctf'], rho, 64, noise,
                        u, np.deg2rad(angle), .5/g['field_A'])
                    signals.append((frame, f'{angle:g}deg_{direction}', signal))
        result.update(noise_std_supplied=noise, field_A=g['field_A'], design_arrays_sha256=sha(ap),
            setup_seconds=time.perf_counter()-start, quadrature_kernel_error=gram.kernel_error,
            preconditioner_diagnostics=gram.preconditioner_diagnostics)
        save()
        for feature in locked['features']:
            name = feature['name']; centers = [feature['center_fraction_field']]
            width = locked['width_fraction_field']; pilot_target = feature['pilot_expected_feature']
            truth = {frame: float(cell_target_coefficients(64, centers, [1], width)@rho)
                for frame, rho in references.items()}
            audit = read(BASE/f'pilot-selected-fixed-v1/{ds}-{name}.json')['targets'][0]['fit']
            for tau in [2., 1.]:
                for pose_model, covariance in [('fixed_pose', None), ('local_gaussian_pose', nuisance)]:
                    tick = time.perf_counter()
                    fit = continuous_gaussian(gram, centers, [1], width, tau, alpha=.05/12,
                        density_radius=2., nuisance=covariance)
                    w = fit.pop('weights'); sd = float(np.linalg.norm(w)); h = fit['half_width']
                    wp = out/f'{ds}-{name}-tau{tau:g}-{pose_model}.npz'
                    np.savez_compressed(wp, weights=w, indices=g['indices'], noise_std=noise)
                    row = dict(target=name, centers=centers, width_fraction_field=width, width_A=20.,
                        pose_model=pose_model, prior_directional_sd=tau, fit=fit,
                        no_data_half_width=2*fit['target_norm'],
                        original_fixed_audit_half_width=audit['half_width'],
                        half_width_over_fixed_audit=h/audit['half_width'],
                        weights_sha256=sha(wp), checks=[])
                    for frame, scenario, signal in signals:
                        mean = float(pilot_target+w@(signal-nominal)); bias = mean-truth[frame]
                        coverage = norm.cdf((h-bias)/sd)-norm.cdf((-h-bias)/sd)
                        power = norm.cdf((np.sign(truth[frame])*mean-h)/sd)
                        row['checks'].append(dict(frame=frame, scenario=scenario, true_target=truth[frame],
                            expected_center=mean, bias=bias, measurement_noise_sd=sd, half_width=h,
                            half_width_over_no_data=h/row['no_data_half_width'],
                            half_width_over_abs_mean=h/abs(mean) if mean else None,
                            analytic_fixed_signal_coverage=float(coverage), correct_sign_probability=float(power)))
                    row['seconds'] = time.perf_counter()-tick
                    result['records'].append(row); save()
                    print(ds, name, tau, pose_model, 'width/audit', row['half_width_over_fixed_audit'],
                        'converged', fit['converged'], 'identity', fit['variance_identity_relative_error'], flush=True)
        result.update(complete=True, seconds=time.perf_counter()-start,
            all_solvers_converged=all(r['fit']['converged'] for r in result['records']),
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
        save()
    except Exception as exc:
        result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--datasets', default='10028,10049,10076')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--output', default='continuous-gaussian-review2-v2')
    args = p.parse_args()
    if args.threads != 1:
        raise ValueError('Version 2 prescribes one CPU thread')
    import finufft
    if args.threads < 1 or (args.threads > 1 and ('tmp/nufft-openmp' not in finufft.__file__ or 'torch' in sys.modules)):
        raise RuntimeError('Use isolated CPU-only threaded FINUFFT')
    files = [Path(__file__), PROTOCOL, ROOT/'src/fourier_splats/uq_continuous_gaussian.py',
        ROOT/'tests/test_continuous_gaussian.py', ROOT/'src/fourier_splats/uq_cached_quadrature.py',
        ROOT/'tests/test_cached_quadrature.py', ROOT/'research/uncertainty/CONTINUOUS-GAUSSIAN-PROTOCOL.md']
    for path in files:
        if subprocess.check_output(['git', 'show', f'HEAD:{path.relative_to(ROOT)}'], cwd=ROOT) != path.read_bytes():
            raise ValueError('Commit protocol and new code before running outcomes')
    snapshot = source_snapshot(ROOT, Path(__file__), [str(p.relative_to(ROOT)) for p in files[1:]]+
        ['scripts/audit_uq_grid_refinement.py', str(LOCK.relative_to(ROOT))])
    rows = read_target_lock(LOCK)['datasets']
    requested = args.datasets.split(',')
    if set(requested)-{r['dataset'] for r in rows}:
        raise ValueError('Dataset outside target lock')
    for locked in rows:
        if locked['dataset'] in requested:
            run(locked, args, snapshot)


if __name__ == '__main__':
    main()
