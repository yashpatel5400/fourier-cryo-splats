#!/usr/bin/env python3
"""Independent 60-decimal Sobolev spot checks; not a global rounding proof."""
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import numpy as np
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_ball_remainder import particle_ball_remainder
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]


def main():
    path = ROOT/'results/uncertainty/development/audit-regressions/enclosure-precision.json'
    if path.exists():
        raise RuntimeError('Preserve existing check')
    source = ROOT/'results/uncertainty/development/continuous-high-band-1024-10A/10049.json'
    fit = json.loads(source.read_text()); cfg = fit['config']; target = fit['targets'][0]
    wp = source.with_name(f"10049-center-{target['width_fraction_field']}-weights.npz")
    saved = np.load(wp)
    g = particle_geometry(ROOT, '10049', 'inference_half0', radius=cfg['frequency_radius'], count=cfg['particles'], seed=cfg['seed'])
    np.testing.assert_array_equal(saved['indices'], g['indices'])
    n, nq = g['ctf'].shape; w = saved['weights'].reshape(n, 2*nq)
    coefficients = (w[:, :nq]+1j*w[:, nq:])*g['ctf']/float(saved['noise_std'])
    mp.mp.dps = 60; rows = []; start = time.perf_counter()
    result = {'complete': False, 'scope': 'Selected-particle high-precision spot checks only, not a uniform numerical certificate',
              'particle_positions': [0, 512, 1023], 'decimal_digits': 60,
              'weights_sha256': hashlib.sha256(wp.read_bytes()).hexdigest(),
              'source_snapshot': source_snapshot(ROOT, Path(__file__)), 'rows': rows}
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    save()
    try:
        for particle in result['particle_positions']:
            k = [[mp.mpf(float(x)) for x in row] for row in g['k'][particle]]
            c = [mp.mpc(float(z.real), float(z.imag)) for z in coefficients[particle]]
            records = {domain: particle_ball_remainder(g['k'][particle], g['q'][particle], coefficients[particle], np.deg2rad(1.), .5/g['field_A'], domain=domain) for domain in ['ball', 'cube']}
            radius = {d: mp.mpf(records[d]['domain_radius']) for d in records}
            terms = {d: [] for d in records}
            def kernel(v, domain):
                r = radius[domain]
                if domain == 'cube':
                    return mp.fprod(2*r if x == 0 else mp.sin(2*mp.pi*r*x)/(mp.pi*x) for x in v)
                z = 2*mp.pi*r*mp.sqrt(mp.fsum(x*x for x in v))
                return 4*mp.pi*r**3/3 if z == 0 else 4*mp.pi*r**3*(mp.sin(z)-z*mp.cos(z))/z**3
            for q in range(nq):
                for l in range(q+1):
                    dot = mp.fsum(x*y for x, y in zip(k[q], k[l]))
                    minus = [x-y for x, y in zip(k[q], k[l])]; plus = [x+y for x, y in zip(k[q], k[l])]
                    for domain in records:
                        value = dot**3*(c[q]*mp.conj(c[l])*kernel(minus, domain)-c[q]*c[l]*kernel(plus, domain))
                        terms[domain].append((1 if q == l else 2)*mp.re(value))
            for domain in records:
                exact = mp.mpf('.5')*(2*mp.pi)**6*mp.fsum(terms[domain])
                diag = records[domain]['sobolev_diagnostics'][2]
                delta = mp.mpf(diag['squared_unpadded'])-exact
                row = {'particle': particle, 'domain': domain, 'order': 3,
                       'high_precision_squared': mp.nstr(exact, 62),
                       'absolute_difference': float(abs(delta)), 'relative_difference': float(abs(delta)/abs(exact)),
                       'heuristic_pad': diag['roundoff_pad'], 'absolute_difference_within_pad': bool(abs(delta) <= diag['roundoff_pad']),
                       'padded_value_above_reference': bool(delta+diag['roundoff_pad'] >= 0)}
                rows.append(row); save(); print(row, flush=True)
        passed = all(r['absolute_difference_within_pad'] and r['padded_value_above_reference'] for r in rows)
        result.update(complete=passed, numerical_failure=not passed, seconds=time.perf_counter()-start); save()
        if not passed:
            raise AssertionError('Precision check failed; outcome retained')
    except Exception as exc:
        result.update(complete=False, error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
