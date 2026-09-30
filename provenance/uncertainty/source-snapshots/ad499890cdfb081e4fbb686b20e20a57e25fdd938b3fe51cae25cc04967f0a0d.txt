#!/usr/bin/env python3
"""Same-weight continuous sign/total-norm class sensitivity, with new centers."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from fourier_splats.uq_data import particle_geometry, VoxelReference
from fourier_splats.uq_continuous import cell_target_coefficients, cell_adjoint, cell_forward
from fourier_splats.uq_sign_class import projected_sign_bounds
from fourier_splats.uq_intervals import bias_aware_half_width_stable, reference_interval_summary
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model
from probe_uq_ball_remainder import validate_source_class

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complete(path):
    result = json.loads(path.read_text())
    if not result.get('complete') or result.get('error'):
        raise ValueError(f'Completed source required: {path}')
    return result


def main():
    ap = BASE/'continuous-high-band-pose-1024-10A/10049-center-1.json'; audit = complete(ap)
    ep = BASE/'expanded-cube-remainder-probe-v2/10049-center-1.json'; enclosure = complete(ep)
    fp = ROOT/audit['config']['fit']; fit_record = complete(fp)
    if sha(fp) != audit['source_fit_sha256'] or sha(ap) != enclosure['source_audit_sha256']:
        raise ValueError('Source certificate changed')
    B, P, _ = validate_source_class(audit, fit_record)
    width = audit['width_fraction_field']; row = next(r for r in fit_record['targets'] if r['target'] == 'center')
    fit = row['fit']; wp = fp.with_name(f'10049-center-{width}-weights.npz')
    if sha(wp) != audit['source_weights_sha256'] or sha(wp) != enclosure['source_weights_sha256']:
        raise ValueError('Source weights changed')
    out = BASE/'sign-class-probe'; out.mkdir(exist_ok=True); path = out/'10049-center-1.json'
    if path.exists():
        raise RuntimeError('Preserve previous outcome')
    start = time.perf_counter()
    result = {'complete': False, 'dataset': '10049', 'target': 'center', 'sigma_A': 10.,
        'source_audit_sha256': sha(ap), 'source_enclosure_sha256': sha(ep), 'source_fit_sha256': sha(fp),
        'source_weights_sha256': sha(wp),
        'scope': 'Alternative origin-norm/negative-part classes, same weights; not a calibrated experimental prior or optimized sign-constrained method',
        'numerical_scope': 'Real-arithmetic support identities with heuristic projection pad; not validated interval arithmetic',
        'source_snapshot': source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py',
            'scripts/probe_uq_ball_remainder.py', 'research/uncertainty/SIGN-CONSTRAINED-DENSITY-PROPOSAL.md']),
        'records': []}
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    save()
    try:
        cfg = fit_record['config']
        g = particle_geometry(ROOT, '10049', 'inference_half0', radius=cfg['frequency_radius'], count=cfg['particles'], seed=cfg['seed'])
        saved = np.load(wp); np.testing.assert_array_equal(saved['indices'], g['indices'])
        w = saved['weights']; noise = float(saved['noise_std'])
        np.testing.assert_allclose(np.linalg.norm(w), audit['noise_sd'], rtol=1e-12)
        cp = BASE/'representation/10049/real_particles-spacing-2.0.npz'
        op, pilot, _, _ = model(g, np.load(cp), 24, noise=noise); pilot = op.expand(pilot)
        pilot_target = float(cell_target_coefficients(24, [[0, 0, 0]], [1], width)@pilot)
        pilot_signal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
        rp = ROOT/'data/uncertainty/references/emd_6487.map'
        ref = VoxelReference.from_mrc(rp, box=64).volume.ravel(); ref /= np.linalg.norm(ref)
        truth = float(cell_target_coefficients(64, [[0, 0, 0]], [1], width)@ref)
        exact_mean = float(w@cell_forward(g['k'], g['ctf'], ref, 64, noise))
        means = {r['scenario']: r['raw_expected_center']-pilot_target+float(w@pilot_signal) for r in audit['reference_checks']}
        np.testing.assert_allclose(means['nominal'], exact_mean, rtol=1e-12, atol=1e-12)
        H = audit['density_bias']/B
        nuisance_field = (audit['pose_polynomial_bias']+enclosure['selected_cubic_bias'])/(B+P)
        negative_norm = float(np.linalg.norm(np.minimum(ref, 0)))
        result.update(reference_sha256=sha(rp), reference_norm=float(np.linalg.norm(ref)),
            reference_negative_part_norm=negative_norm, reference_target=truth,
            continuous_residual_norm_upper=H, nuisance_field_norm_upper=nuisance_field,
            noise_sd=audit['noise_sd'], alpha_noise=audit['alpha_noise'], delta_spectral=audit['alpha_numerical'],
            interval_scope='Pointwise per declared class; not a guarantee after selecting a class using these results.')
        for box in [24, 48, 96]:
            ell = cell_target_coefficients(box, [[0, 0, 0]], [1], width)
            adjoint = cell_adjoint(g['k'], g['ctf'], w, box, noise)
            projection = ell-adjoint
            error = 1e-8*(1+np.linalg.norm(ell)+np.linalg.norm(adjoint))
            pp = out/f'10049-center-projection-{box}.npz'
            np.savez(pp, coefficients=projection)
            for radius in [1., 2., 3.]:
                for fraction in [0., .1, .25, .5, 1.]:
                    eta = fraction*radius
                    bound = projected_sign_bounds(projection, H, radius, eta, error)
                    no_center = (radius-eta)*fit['target_norm']/2
                    no_half = (radius+eta)*fit['target_norm']/2
                    item = {'box': box, 'radius': radius, 'negative_radius': eta,
                        'source_projection_sha256': sha(pp), 'support': bound,
                        'reference_in_class': bool(np.linalg.norm(ref) <= radius+1e-12 and negative_norm <= eta),
                        'no_data_center': no_center, 'no_data_half_width': no_half, 'checks': []}
                    for pose_label, field_bound, scenarios in [('fixed_pose', 0., ['nominal']),
                            ('one_degree_half_A', nuisance_field, list(means))]:
                        bias = bound['bias_half_width']+radius*field_bound
                        half = bias_aware_half_width_stable(audit['noise_sd'], bias, audit['alpha_noise'])
                        for scenario in scenarios:
                            expected = means[scenario]+bound['center_offset']
                            check = reference_interval_summary(truth, expected, audit['noise_sd'], no_center, half, no_half)
                            item['checks'].append(dict(pose_class=pose_label, scenario=scenario,
                                bias_upper=bias, raw_half_width=half, selected_half_width=min(half, no_half),
                                relative_half_width=min(half, no_half)/no_half, **check))
                    result['records'].append(item); save()
            print('DONE_PROJECTION', box, 'reference negative norm', negative_norm, flush=True)
        result.update(complete=True, seconds=time.perf_counter()-start); save()
    except Exception as exc:
        result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
