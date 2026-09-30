#!/usr/bin/env python3
"""Execute the two prespecified remainder-only feasibility cases, without fitting."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_higher_remainder import higher_order_particle_remainders
from fourier_splats.uq_provenance import source_snapshot
from probe_uq_ball_remainder import validate_source_class

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
CASES = ['pilot-selected-pose-v1/10049-pilot_region_1-1.json',
         'continuous-high-band-pose-1024-10A/10049-center-1.json']
DEGREES = (2, 3, 4, 5)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    out = BASE/'higher-order-remainder-probe'
    out.mkdir(parents=True, exist_ok=True)
    if any(out.glob('*.json')):
        raise RuntimeError('Preserve existing outcomes; no overwriting or outcome-based continuation')
    snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/probe_uq_ball_remainder.py',
        'research/uncertainty/HIGHER-ORDER-REMAINDER-PROBE.md'])
    errors = []
    for name in CASES:
        start = time.perf_counter(); audit_path = BASE/name; path = out/audit_path.name
        result = {'complete': False, 'source_audit': str(audit_path.relative_to(ROOT)),
            'source_snapshot': snapshot, 'degrees': DEGREES,
            'scope': 'Remainder-only feasibility diagnostic; no higher-order interval, coverage or power result',
            'numerical_scope': 'Real-arithmetic identities with heuristic floating-point guards, not interval arithmetic'}
        def save():
            path.write_text(json.dumps(result, indent=2)+'\n')
        save()
        try:
            audit = json.loads(audit_path.read_text())
            if not audit.get('complete') or audit.get('error'):
                raise ValueError('Completed source audit required')
            fit_path = ROOT/audit['config']['fit']; fit = json.loads(fit_path.read_text())
            if sha(fit_path) != audit['source_fit_sha256']:
                raise ValueError('Source fit changed')
            B, P, origin = validate_source_class(audit, fit)
            if (B, P) != (2., 1.) or audit['config']['angle'] != 1. or audit['config']['shift_A'] != .5:
                raise ValueError('The prespecified case uses B=2, P=1, one degree and 0.5 A')
            ds, target, width = audit['dataset'], audit['target'], audit['width_fraction_field']
            wp = fit_path.with_name(f'{ds}-{target}-{width}-weights.npz')
            if sha(wp) != audit['source_weights_sha256']:
                raise ValueError('Source weights changed')
            cfg = fit['config']
            g = particle_geometry(ROOT, ds, 'inference_half0', radius=cfg['frequency_radius'], count=cfg['particles'], seed=cfg['seed'])
            saved = np.load(wp); np.testing.assert_array_equal(saved['indices'], g['indices'])
            n, nq = g['ctf'].shape; weights = saved['weights'].reshape(n, 2*nq)
            coefficients = (weights[:, :nq]+1j*weights[:, nq:])*g['ctf']/float(saved['noise_std'])
            records = []
            for i in range(n):
                records.append(higher_order_particle_remainders(g['k'][i], g['q'][i], coefficients[i],
                    np.deg2rad(1.), .5/g['field_A'], degrees=DEGREES, domain='cube'))
            field = np.array([[r['field_remainder'] for r in row['records']] for row in records])
            residual = np.array([[r['embedding_residual_remainder'] for r in row['records']] for row in records])
            np.savez(path.with_suffix('.npz'), taylor_degrees=DEGREES, field_remainders=field,
                embedding_residual_remainders=residual, sobolev_norms=np.array([r['sobolev_norms'] for r in records]), indices=g['indices'])
            totals = (B+P)*field.sum(axis=0)
            original = (B+P)*sum(r['degree_two_existing_bound'] for r in records)
            np.testing.assert_allclose(totals[0], original, rtol=5e-13, atol=1e-12)
            result.update(complete=True, dataset=ds, target=target, particles=n, width_A=width*g['field_A'],
                source_audit_sha256=sha(audit_path), source_fit_sha256=sha(fit_path), source_weights_sha256=sha(wp),
                density_radius=B, pilot_norm_bound=P, class_metadata_origin=origin,
                domain='cube', joint_ball=True, rotation_radius_degrees=1., translation_radius_A=.5,
                remainder_bias_components=[{'taylor_degree': degree, 'remainder_bias_upper': float(value),
                    'relative_to_degree_two': float(value/totals[0]), 'higher_pose_fields_implemented': degree == 2}
                    for degree, value in zip(DEGREES, totals)],
                maximum_degree_two_absolute_difference=max(r['degree_two_absolute_difference'] for r in records),
                maximum_relative_sobolev_pad=max(r['maximum_relative_sobolev_pad'] for r in records),
                embedding_residual_bias_components=((B+P)*residual.sum(axis=0)).tolist(),
                arrays_sha256=sha(path.with_suffix('.npz')), seconds=time.perf_counter()-start)
            save(); print('COMPLETE', name, totals.tolist(), 'seconds', result['seconds'], flush=True)
        except Exception as exc:
            result.update(error=repr(exc), seconds=time.perf_counter()-start); save()
            errors.append(repr(exc)); print('FAILED_RETAINED', name, repr(exc), flush=True)
    if errors:
        raise RuntimeError(f'Incomplete diagnostic; both declared cases attempted: {errors}')


if __name__ == '__main__':
    main()
