#!/usr/bin/env python3
"""Frozen paired local-alignment simulation; every failure and raw interval saved."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import norm
from fourier_splats.uq_continuous import cell_forward, cell_target_coefficients, continuous_certificate
from fourier_splats.uq_continuous_gaussian import continuous_gaussian, cell_pose_jacobian
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_fourier_pose_baseline import GaussianNuisanceWhitening
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_local_alignment import CachedCellTemplate, refine_local_poses
from fourier_splats.uq_mixed_pose import mixed_pose_interval
from fourier_splats.uq_end_to_end import dephase_observations, observed_interval
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
PROTOCOL = ROOT/'research/uncertainty/END-TO-END-LOCAL-POSE-PROTOCOL.md'
CALIBRATION = BASE/'local-alignment-calibration-v1'
TARGETS = [('center', [[0., 0., 0.]], [1.]),
           ('contrast', [[0., 0., .08], [0., 0., -.08]], [1., -1.])]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def calibration_parameters(ds):
    directory = CALIBRATION/ds
    summary = json.loads((directory/'summary.json').read_text())
    if not summary.get('complete') or summary.get('error') or len(summary['records']) != 128:
        raise ValueError('Completed 128-dataset calibration required before any study coverage')
    hashes = {'summary.json': sha(directory/'summary.json')}
    parameters = {}
    for template in ['oracle_reference', 'independent_pilot']:
        errors, scores = [], []
        for entry in summary['records']:
            rep = entry['replicate']; rp = directory/f'replicate-{rep:03d}.json'
            if sha(rp) != entry['record_sha256']:
                raise ValueError('Calibration record changed')
            record = json.loads(rp.read_text()); ap = directory/f'replicate-{rep:03d}.npz'
            if not record.get('complete') or sha(ap) != record['arrays_sha256']:
                raise ValueError('Incomplete or changed calibration arrays')
            hashes[rp.name] = sha(rp); hashes[ap.name] = sha(ap)
            with np.load(ap) as saved:
                e = np.concatenate([saved[template+'_rotation_error_vector_degrees'], saved[template+'_shifts_A']], axis=1)
            errors.append(e)
            scores.append(float(np.linalg.norm(e/np.array([1.]*3+[.5]*2), axis=1).max()))
        error = np.stack(errors)
        parameters[template] = dict(joint_score_bound=max(scores), calibration_dataset_scores=scores,
            gaussian_coordinate_rms=np.sqrt(np.mean(error**2, axis=(0, 1))).tolist(),
            mean_error_coordinates=error.mean(axis=(0, 1)).tolist(),
            per_particle_mean_error_coordinates=error.mean(axis=0).tolist(),
            known_truth_tolerance_failure_probability=1/129,
            scope='Same fixed simulated truth/design/optimizer exchangeability only; not experimental pose calibration.')
    gp = directory/'generators.npz'
    if sha(gp) != summary['generators_sha256']:
        raise ValueError('Calibration generator changed')
    hashes[gp.name] = sha(gp)
    return parameters, hashes, gp


def run(ds, args, snapshot):
    out = BASE/args.output/ds
    params, cal_hashes, gp = calibration_parameters(ds)
    if out.exists():
        raise RuntimeError(f'Preserve previous study: {out}')
    out.mkdir(parents=True)
    start = time.perf_counter(); summary_path = out/'summary.json'
    summary = dict(complete=False, dataset=ds, protocol_sha256=sha(PROTOCOL), source_snapshot=snapshot,
        config=dict(replicates=200, particles=128, radius=5, seed_root=261002,
            width_fraction_field=.07, density_radius=2., alpha=.05, pose_audit_noise_alpha=.04,
            quadrature_order=40, preconditioner_rank=1024, threads=args.threads),
        calibration_parameters=params, calibration_hashes=cal_hashes, records=[],
        scope='Paired local-alignment simulation with a fixed generator/template. No experimental coverage or global unknown-pose reconstruction claim.')
    def save():
        summary_path.write_text(json.dumps(summary, indent=2)+'\n')
    save()
    with np.load(gp) as data:
        g = {k: data[k].copy() for k in ['k', 'q', 'ctf', 'indices']}
        rho, pilot = data['truth'].copy(), data['pilot'].copy()
        noise, field = float(data['noise_std']), float(data['field_A'])
        signal = data['noiseless_signal'].copy()
    templates = [('oracle_reference', CachedCellTemplate(rho, 64)),
                 ('independent_pilot', CachedCellTemplate(pilot, 24))]
    truths = {name: float(cell_target_coefficients(64, c, s, .07)@rho) for name, c, s in TARGETS}
    pilot_targets = {name: float(cell_target_coefficients(24, c, s, .07)@pilot) for name, c, s in TARGETS}
    for replicate in range(200):
        tick = time.perf_counter(); rp = out/f'replicate-{replicate:03d}.json'
        row = dict(complete=False, dataset=ds, replicate=replicate, templates=[], intervals=[])
        arrays = {}
        def save_row():
            rp.write_text(json.dumps(row, indent=2)+'\n')
        save_row()
        try:
            rng = np.random.default_rng(np.random.SeedSequence([261002, int(ds), replicate]))
            y_a = signal+rng.normal(size=signal.shape)
            r0 = Rotation.from_rotvec(rng.normal(size=(128, 3))*np.deg2rad(2.)/np.sqrt(3)).as_matrix()
            t0 = rng.normal(size=(128, 2))*.5/np.sqrt(2)
            y_b = signal+rng.normal(size=signal.shape)
            arrays.update(y_a=y_a, y_b=y_b, initial_rotation=r0, initial_shift_A=t0)
            for template_name, template in templates:
                align = refine_local_poses(template, g['k'], g['q'], g['ctf'], y_a, noise, field,
                    initial_rotation=r0, initial_shift_A=t0)
                rotations, shifts = align.pop('rotations'), align.pop('shifts_A')
                rot_error = Rotation.from_matrix(rotations).as_rotvec()*180/np.pi
                observed_score = float(np.linalg.norm(np.concatenate([rot_error, shifts/.5], axis=1), axis=1).max())
                for name, value in [('rotations', rotations), ('shifts_A', shifts),
                        ('initial_objective', align.pop('initial_objective')), ('final_objective', align.pop('final_objective'))]:
                    arrays[template_name+'_'+name] = value
                calibration = params[template_name]; bound = calibration['joint_score_bound']
                angle, shift_fraction = np.deg2rad(bound), .5*bound/field
                alignment_row = dict(template=template_name, alignment=align, observed_joint_score=observed_score,
                    pose_inside_calibrated_ball=bool(observed_score <= bound),
                    rotation_rms_degrees=float(np.sqrt(np.mean(np.sum(rot_error**2, axis=1)))),
                    shift_rms_A=float(np.sqrt(np.mean(np.sum(shifts**2, axis=1)))), fits=[])
                row['templates'].append(alignment_row)
                k = g['k']@rotations
                gram = QuadratureObservationGram(k, g['ctf'], noise, order=40, preconditioner_rank=1024)
                gram.nthreads = args.threads
                prior_prediction = cell_forward(k, g['ctf'], pilot, 24, noise)
                observed = {name: dephase_observations(y, g['q'], shifts, field) for name, y in
                    [('same_image', y_a), ('independent_image', y_b)]}
                jac = cell_pose_jacobian(k, g['q'], g['ctf'], pilot, 24, noise)
                rms = np.asarray(calibration['gaussian_coordinate_rms'])*np.array([np.pi/180]*3+[1/field]*2)
                nuisance = GaussianNuisanceWhitening(jac*rms)
                for target_name, centers, signs in TARGETS:
                    prefix = template_name+'_'+target_name
                    fit = continuous_certificate(gram, centers, signs, .07, 2., alpha=.05, maxiter=100, rtol=.005)
                    w = fit.pop('weights'); arrays[prefix+'_audit_weights'] = w
                    sd = float(np.linalg.norm(w)); residual = fit['bias']/2
                    pair = np.einsum('nmp,nm->np', jac*np.array([angle]*3+[shift_fraction]*2), w.reshape(128, -1))
                    mixed = mixed_pose_interval(k, g['q'], g['ctf'], w, noise, angle, shift_fraction,
                        2., 1., residual, alpha=.05, pilot_derivative_pairings=pair)
                    # The same first-order norm/remainder calculation serves the deterministic class.
                    deterministic_bias = mixed['density_bias']+sum(mixed['scalar_derivative_bounds'])+mixed['nonlinear_remainder_bias']
                    common_half = (mixed['density_bias']+.1*sum(mixed['scalar_derivative_bounds'])
                        +1.1**2*mixed['nonlinear_remainder_bias']+mixed['subgaussian_tail'])
                    methods = [('fixed_folded', w, bias_aware_half_width_stable(sd, fit['bias'], .05), True),
                        ('fixed_sum', w, fit['bias']+norm.isf(.025)*sd, True),
                        ('deterministic_pose', w, bias_aware_half_width_stable(sd, deterministic_bias, .04), True),
                        ('mixed_common0', w, mixed['half_width'], True),
                        ('mixed_common01', w, common_half, True)]
                    fit_record = dict(target=target_name, fixed_audit=fit, mixed_audit=mixed,
                        deterministic_pose_bias=float(deterministic_bias), gaussian_fits=[])
                    alignment_row['fits'].append(fit_record)
                    for tau in [2., 1.]:
                        for label, cov in [('fixed', None), ('pose', nuisance)]:
                            gaussian = continuous_gaussian(gram, centers, signs, .07, tau, alpha=.05,
                                density_radius=2., nuisance=cov)
                            gw = gaussian.pop('weights'); name = f'gaussian_tau{tau:g}_{label}'
                            arrays[prefix+'_'+name+'_weights'] = gw
                            fit_record['gaussian_fits'].append(dict(method=name, **gaussian))
                            methods.append((name, gw, gaussian['half_width'], False))
                    no_data = 2*fit['target_norm']
                    for method, weights, half, fallback in methods:
                        for image_mode, y in observed.items():
                            center = float(pilot_targets[target_name]+weights@(y-prior_prediction))
                            row['intervals'].append(dict(template=template_name, target=target_name,
                                method=method, image_mode=image_mode,
                                **observed_interval(center, half, pilot_targets[target_name], no_data, truths[target_name], fallback)))
                    save_row()
                del gram
            ap = out/f'replicate-{replicate:03d}.npz'
            np.savez_compressed(ap, **arrays)
            row.update(complete=True, arrays_sha256=sha(ap), seconds=time.perf_counter()-tick)
            save_row()
        except Exception as exc:
            # Save the failure in the planned denominator, retain partial arrays, continue.
            ap = out/f'replicate-{replicate:03d}-partial.npz'
            np.savez_compressed(ap, **arrays)
            row.update(error=repr(exc), partial_arrays_sha256=sha(ap), seconds=time.perf_counter()-tick)
            save_row()
        summary['records'].append(dict(replicate=replicate, record_sha256=sha(rp),
            complete=row['complete'], error=row.get('error'), seconds=row['seconds']))
        save()
        print(ds, replicate, 'complete', row['complete'], 'seconds', row['seconds'], 'error', row.get('error'), flush=True)
    summary.update(complete=True, all_replicates_successful=all(r['complete'] for r in summary['records']),
        seconds=time.perf_counter()-start,
        peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
    save()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--datasets', default='10028,10049,10076')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--output', default='end-to-end-local-pose-v1')
    args = parser.parse_args()
    if args.threads != 1:
        raise ValueError('Frozen initial study uses one NUFFT thread per process')
    files = [Path(__file__), PROTOCOL, ROOT/'src/fourier_splats/uq_end_to_end.py', ROOT/'tests/test_end_to_end.py',
        ROOT/'src/fourier_splats/uq_mixed_pose.py', ROOT/'src/fourier_splats/uq_local_alignment.py',
        ROOT/'src/fourier_splats/uq_continuous_gaussian.py']
    for path in files:
        if subprocess.check_output(['git', 'show', f'HEAD:{path.relative_to(ROOT)}'], cwd=ROOT) != path.read_bytes():
            raise ValueError('Commit protocol and sources before study outcomes')
    datasets = args.datasets.split(',')
    if set(datasets)-{'10028', '10049', '10076'}:
        raise ValueError('Undeclared geometry')
    # Check all requested calibrations first, so a partially ready batch cannot start.
    for ds in datasets:
        if not json.loads((CALIBRATION/ds/'summary.json').read_text()).get('complete'):
            raise ValueError('Calibration not complete')
    snapshot = source_snapshot(ROOT, Path(__file__), [str(p.relative_to(ROOT)) for p in files[1:]])
    for ds in datasets:
        run(ds, args, snapshot)


if __name__ == '__main__':
    main()
