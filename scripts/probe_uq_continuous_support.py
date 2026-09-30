#!/usr/bin/env python3
"""Declared-support development at a 10 Angstrom Gaussian target scale.

This changes the assumed class, does not estimate it, and retains full-map
outside energy. All reference outcomes are diagnostics, never fit criteria.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry, VoxelReference
from fourier_splats.uq_continuous import continuous_certificate, cell_forward, cell_target_coefficients
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_continuous_pose import pose_cell_forward
from fourier_splats.uq_continuous_moments import integrated_cubic_remainder
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_random_spectral import gaussian_power_upper
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_support import rescale_supported_problem, restrict_cell_density, outside_tail_bias_upper
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'
MAPS = {'10028': '2660', '10049': '6487', '10076': '8434'}


def run(args, dataset, snapshot):
    out = BASE/args.output; out.mkdir(exist_ok=True)
    g = particle_geometry(ROOT, dataset, 'inference_half0', radius=args.radius,
                          count=args.particles, seed=609315)
    noise = float(np.load(BASE/'continuous-quadrature-optimized'/f'{dataset}-center-0.07-weights.npz')['noise_std'])
    ck = np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
    op, pilot, _, _ = model(g, ck, 24, noise=noise); pilot = op.expand(pilot)
    rho = VoxelReference.from_mrc(ROOT/f'data/uncertainty/references/emd_{MAPS[dataset]}.map', box=64).volume.ravel()
    rho /= np.linalg.norm(rho)
    full_gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=80, preconditioner_rank=0)
    width = args.width_A/g['field_A']; B = 2.; angle = np.deg2rad(args.angle)
    shift = args.shift_A/g['field_A']; delta = 1e-6
    a_full, norm2_full = full_gram.target([[0, 0, 0]], [1], width)
    for side in map(float, args.sides.split(',')):
        path = out/f'{dataset}-side{side:g}.json'
        if path.exists(): raise RuntimeError('Preserve existing support outcomes')
        start = time.perf_counter()
        p = rescale_supported_problem(g['k'], g['q'], g['ctf'], [[0, 0, 0]], [1], width, shift, side)
        pc, pb, pilot_tail = restrict_cell_density(pilot, 24, side)
        rc, rb, ref_tail = restrict_cell_density(rho, 64, side)
        P = float(np.linalg.norm(pc))
        record = {'stage': 'Exploratory support ablation; known simulated noise and supplied pose/density/support bounds',
                  'config': vars(args), 'dataset': dataset, 'side_fraction_field': side,
                  'support_side_A': side*g['field_A'], 'source_snapshot': snapshot,
                  'noise_std_fixed_to_radius5_design': noise, 'target_sigma_A': args.width_A,
                  'density_radius_inside': B, 'pilot_inside_norm': P,
                  'original_pilot_outside_norm': pilot_tail, 'full_reference_outside_norm': ref_tail,
                  'complete': False, 'fit_history': [], 'spectral_history': [], 'records': []}
        def save(): path.write_text(json.dumps(record, indent=2)+'\n')
        save()
        try:
            order = int(np.ceil(32+48*side))
            gram = QuadratureObservationGram(p['k'], p['ctf'], noise, order=order,
                                             preconditioner_rank=args.preconditioner_rank)
            def progress(row):
                record['fit_history'].append(row); save()
                print(dataset, side, 'FIT', row['iteration'], row['relative_gap'], flush=True)
            fit = continuous_certificate(gram, p['centers'], p['signs'], p['width'], B,
                                         maxiter=args.maxiter, rtol=.005, callback=progress)
            w = fit.pop('weights'); fit.pop('history'); record['fit'] = fit
            np.savez(out/f'{dataset}-side{side:g}-weights.npz', weights=w, indices=g['indices'], noise_std=noise)
            del gram
            # The fixed-pose fitting weights precede fresh audit probes.
            field = PolynomialPoseFieldOperator(p['k'], p['q'], p['ctf'], w, noise, angle, p['shift'], order=order)
            scaling = field.establish_coefficient_scaling()
            def spectral_progress(row):
                record['spectral_history'].append(row); save()
                print(dataset, side, 'SPECTRAL', row['iteration'], row['upper_over_lower'], flush=True)
            spectral = gaussian_power_upper(field.spatial_gram, field.shape[0], delta, 4, args.power_iterations,
                609801+int(dataset)*100+int(side*100), scaling['quadrature_gram_trace_upper'], callback=spectral_progress)
            pose_field = min(scaling['triangle_field_bound'], np.sqrt(scaling['sum_group_scales']*(
                spectral['eigenvalue_upper']+scaling['integration_eigenvalue_pad'])))
            del field
            moment = integrated_cubic_remainder(p['k'], p['q'], p['ctf'], w, noise, angle, p['shift'], B, P)
            inside_bias = fit['bias']+(B+P)*pose_field+moment['remainder_bias']
            gw = full_gram.matvec(w); errors = full_gram.quadrature_error(w)
            roundoff_pad = 50*np.finfo(float).eps*len(w)*(norm2_full+2*abs(w@a_full)+abs(w@gw))
            full_residual = np.sqrt(max(0., norm2_full-2*w@a_full+w@gw+errors['squared_field_norm'])+roundoff_pad)
            outside = outside_tail_bias_upper(g['k'], g['q'], g['ctf'], w, noise, angle, shift, full_residual, 1.)
            record.update(pose_field_bound=float(pose_field), inside_bias=float(inside_bias),
                          cubic_bias=moment['remainder_bias'], unit_tail_bias=outside,
                          spectral_certificate=spectral, alpha_noise=.05-delta, alpha_numerical=delta)
            nominal = cell_forward(p['k'], p['ctf'], pc, pb, noise)
            pilot_target = float(cell_target_coefficients(pb, p['centers'], p['signs'], p['width'])@pc)
            rng = np.random.default_rng(609811+int(dataset))
            signals = []
            for kind in ['cropped_reference', 'full_reference']:
                actual_g = p if kind == 'cropped_reference' else g
                actual_rho, box = (rc, rb) if kind == 'cropped_reference' else (rho, 64)
                actual_shift = p['shift'] if kind == 'cropped_reference' else shift
                actual_centers, actual_signs, actual_width = ((p['centers'], p['signs'], p['width'])
                    if kind == 'cropped_reference' else ([[0, 0, 0]], [1], width))
                truth = float(cell_target_coefficients(box, actual_centers, actual_signs, actual_width)@actual_rho)
                for scenario in ['nominal', 'coherent_x', 'random_boundary']:
                    u = np.zeros((len(g['k']), 5))
                    if scenario == 'coherent_x': u[:, 0] = 1.
                    if scenario == 'random_boundary':
                        u = rng.normal(size=u.shape); u /= np.linalg.norm(u, axis=1)[:, None]
                    signal = pose_cell_forward(actual_g['k'], actual_g['q'], actual_g['ctf'], actual_rho,
                                               box, noise, u, angle, actual_shift)
                    signals.append((kind, scenario, truth, float(pilot_target+w@(signal-nominal))))
            for tail in map(float, args.tail_radii.split(',')):
                bias = inside_bias+tail*outside['bias_upper']
                raw_half = bias_aware_half_width_stable(np.linalg.norm(w), bias, .05-delta)
                no_data = B*fit['target_norm']+tail*np.sqrt(norm2_full)
                fallback = raw_half >= no_data; half = min(raw_half, no_data)
                row = {'tail_radius_supplied': tail, 'half_width_before_fallback': float(raw_half),
                       'relative_half_width': float(half/no_data), 'uses_no_data': bool(fallback),
                       'reference_checks': []}
                for kind, scenario, truth, raw_mean in signals:
                    mean = pilot_target if fallback else raw_mean; sd = 0. if fallback else float(np.linalg.norm(w))
                    actual_bias = mean-truth
                    coverage = (float(norm.cdf((half-actual_bias)/sd)-norm.cdf((-half-actual_bias)/sd))
                                if sd else float(abs(actual_bias) <= half))
                    power = float(norm.cdf((np.sign(truth)*mean-half)/sd)) if sd else float(np.sign(truth)*mean > half)
                    row['reference_checks'].append({'kind': kind, 'scenario': scenario, 'target': truth,
                        'expected_center': mean, 'analytic_coverage': coverage, 'correct_sign_probability': power,
                        'sufficient_class_membership_check': bool(kind == 'cropped_reference' or ref_tail <= tail)})
                record['records'].append(row); save()
                print(dataset, side, 'TAIL', tail, 'WIDTH', row['relative_half_width'], flush=True)
            record.update(complete=True, seconds=time.perf_counter()-start); save()
        except Exception as exc:
            record.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--datasets', default='10049'); p.add_argument('--sides', default='.5,.75,1')
    p.add_argument('--particles', type=int, default=128); p.add_argument('--radius', type=float, default=12)
    p.add_argument('--width-A', type=float, default=10); p.add_argument('--angle', type=float, default=1)
    p.add_argument('--shift-A', type=float, default=.5); p.add_argument('--tail-radii', default='0,.01,.05,.1')
    p.add_argument('--preconditioner-rank', type=int, default=2048); p.add_argument('--maxiter', type=int, default=100)
    p.add_argument('--power-iterations', type=int, default=24); p.add_argument('--output', default='continuous-support-probe')
    args = p.parse_args(); snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py'])
    for dataset in args.datasets.split(','): run(args, dataset, snapshot)


if __name__ == '__main__': main()
