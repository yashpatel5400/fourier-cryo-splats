#!/usr/bin/env python3
"""One numerical-only radius-12 case; no image or coverage outcomes."""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from scipy.sparse.linalg import cg, LinearOperator
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_gaussian_preconditioner import corrected_gram_preconditioner
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_provenance import source_snapshot


def main():
    root = Path(__file__).resolve().parents[1]
    out = root/'results/uncertainty/development/audit-regressions/continuous-gaussian-diagonal-preflight.json'
    if out.exists():
        raise RuntimeError('Preserve previous preflight')
    record = dict(scope='Numerical residuals only; same first 10028 locked target, continuous operator and tau=2. No pixels, reference density outcomes or coverage. Standard residual-diagonal low-rank preconditioner; old preconditioner and failed solves preserved.',
        source_snapshot=source_snapshot(root, Path(__file__), ['src/fourier_splats/uq_gaussian_preconditioner.py', 'tests/test_gaussian_preconditioner.py']), records=[])
    d = json.loads((root/'results/uncertainty/development/pilot-selected-fixed-v1/10028-pilot_region_1.json').read_text())
    t = d['targets'][0]
    wp = root/f'results/uncertainty/development/pilot-selected-fixed-v1/10028-pilot_region_1-{t["width_fraction_field"]}-weights.npz'
    digest = hashlib.sha256(wp.read_bytes()).hexdigest()
    if digest != t['source_weights_sha256']:
        raise ValueError('Archived initialization changed')
    data = np.load(wp); noise = float(data['noise_std']); w0 = data['weights']
    g = particle_geometry(root, '10028', 'inference_half0', radius=12, count=128, seed=609315)
    tick = time.perf_counter()
    gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=80, preconditioner_rank=2048)
    gram.nthreads = 1
    pre = corrected_gram_preconditioner(gram, .25)
    record.update(setup_seconds=time.perf_counter()-tick, initialization_sha256=digest,
        gram_preconditioner_diagnostics=gram.preconditioner_diagnostics)
    a, _ = gram.target(t['centers_fraction_field'], [1], t['width_fraction_field'])
    system = LinearOperator(gram.shape, matvec=lambda x: gram.matvec(x)+.25*x, dtype=float)
    for name, initial in [('zero', np.zeros_like(a)), ('archived_audit', w0)]:
        count = [0]; trace = []; start = time.perf_counter()
        def callback(w):
            count[0] += 1
            if count[0] % 25 == 0:
                residual = float(np.linalg.norm(system@w-a)/np.linalg.norm(a))
                trace.append([count[0], residual]); print(name, count[0], residual, flush=True)
        w, info = cg(system, a, x0=initial, M=pre, rtol=1e-10, atol=0., maxiter=300, callback=callback)
        record['records'].append(dict(initialization=name, info=int(info), iterations=count[0],
            seconds=time.perf_counter()-start, relative_residual=float(np.linalg.norm(system@w-a)/np.linalg.norm(a)), trace=trace))
        out.write_text(json.dumps(record, indent=2)+'\n')


if __name__ == '__main__':
    main()
