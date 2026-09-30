#!/usr/bin/env python3
"""One compute-only probe of a coordinate metric; no weight optimization."""
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_sobolev_penalty import SobolevRemainderPenalty
from fourier_splats.uq_block_preconditioner import cubic_remainder_metric
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    source = BASE/'pilot-selected-fixed-v1/10049-pilot_region_1.json'
    fixed = json.loads(source.read_text()); row = fixed['targets'][0]
    if not fixed.get('complete') or fixed.get('error'):
        raise ValueError('Completed fixed-pose source required')
    wp = source.with_name(f'10049-pilot_region_1-{row["width_fraction_field"]}-weights.npz')
    out = BASE/'cubic-preconditioner-probe'; out.mkdir(exist_ok=True)
    path = out/'10049-pilot_region_1.json'
    if path.exists():
        raise RuntimeError('Preserve previous diagnostic')
    result = {'complete': False, 'scope': 'Compute-only coordinate metric at original fixed-pose weights; no new estimator, spectral certificate or confidence interval.',
        'source_fit_sha256': sha(source), 'source_weights_sha256': sha(wp),
        'config': {'rank': 1024, 'seed': 650171, 'probes': 4, 'relative_floor': 1e-8,
                   'relative_check_threshold': 1e-5, 'angle_degrees': 1., 'shift_A': .5},
        'source_snapshot': source_snapshot(ROOT, Path(__file__), ['research/uncertainty/CUBIC-PRECONDITIONER-DEVELOPMENT.md'])}
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    save(); start = time.perf_counter()
    try:
        saved = np.load(wp); w = saved['weights']; noise = float(saved['noise_std'])
        g = particle_geometry(ROOT, '10049', 'inference_half0', radius=12, count=128, seed=fixed['config']['seed'])
        np.testing.assert_array_equal(saved['indices'], g['indices'])
        gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=80, preconditioner_rank=1024)
        penalty = SobolevRemainderPenalty(g['k'], g['q'], g['ctf']/noise, np.deg2rad(1.), .5/g['field_A'])
        result['input_setup_seconds'] = time.perf_counter()-start; stamp = time.perf_counter()
        H = row['fit']['bias']/2; z = float(norm.isf((.05/12-1e-6/12)/2))
        transform = cubic_remainder_metric(gram, penalty, w, 2., 1., H, z)
        result.update(metric_setup_seconds=time.perf_counter()-stamp, diagnostics=transform.diagnostics)
        save(); print('METRIC_READY', result['metric_setup_seconds'], flush=True)
        # Independently apply the stated metric, not the factors used for roots.
        ww = w.reshape(128, 2, penalty.nq); norms = np.empty((128, 4))
        for i in range(128):
            for d in range(4):
                norms[i, d] = np.sqrt(max(0., sum(ww[i, c]@penalty.blocks[i, d, c]@ww[i, c] for c in range(2))))
        floor = np.maximum(1e-8*np.median(norms, axis=0), 1e-30)
        multipliers = 3*penalty.multipliers/np.maximum(norms, floor)
        def metric_action(v):
            packed = v.reshape(128, 2, penalty.nq); out = (z/np.linalg.norm(w))*packed
            for i in range(128):
                for d in range(4):
                    for c in range(2):
                        out[i, c] += multipliers[i, d]*(penalty.blocks[i, d, c]@packed[i, c])
            for c, (factor, _) in enumerate(gram.preconditioner_factors):
                out[:, c] += ((2/H)*factor@(factor.T@packed[:, c].ravel())).reshape(128, penalty.nq)
            return out.ravel()
        rng = np.random.default_rng(650171); checks = []; stamp = time.perf_counter()
        for j in range(4):
            x = rng.normal(size=w.size); x /= np.linalg.norm(x)
            y = rng.normal(size=w.size); y /= np.linalg.norm(y)
            tx = transform.to_weights(x)
            checks.append({'probe': j, 'roundtrip_relative_error': float(np.linalg.norm(transform.to_coordinates(tx)-x)),
                'whitening_relative_error': float(np.linalg.norm(transform.gradient_to_coordinates(metric_action(tx))-x)),
                'gradient_duality_absolute_error': float(abs(y@tx-transform.gradient_to_coordinates(y)@x))})
        result.update(complete=True, checks=checks, check_seconds=time.perf_counter()-stamp,
            checks_passed=all(r['roundtrip_relative_error'] < 1e-5 and r['whitening_relative_error'] < 1e-5 and r['gradient_duality_absolute_error'] < 1e-8 for r in checks),
            seconds=time.perf_counter()-start,
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
        save(); print('COMPLETE', result['checks_passed'], result['seconds'], checks, flush=True)
    except Exception as exc:
        result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
