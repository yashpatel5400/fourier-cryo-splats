#!/usr/bin/env python3
"""One high-rank numerical preflight, no new pixels/reference coverage."""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from scipy.sparse.linalg import cg, LinearOperator
from fourier_splats.uq_cached_quadrature import CachedQuadratureGram
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_provenance import source_snapshot


def main():
    root = Path(__file__).resolve().parents[1]; base = root/'results/uncertainty/development'
    out = base/'audit-regressions/continuous-gaussian-rank8192-preflight.json'
    if out.exists():
        raise RuntimeError('Preserve previous probe')
    source = json.loads((base/'continuous-gaussian-review2-v1/10028.json').read_text())['records'][0]
    wp = base/'continuous-gaussian-review2-v1/10028-pilot_region_1-tau2-fixed_pose.npz'
    if hashlib.sha256(wp.read_bytes()).hexdigest() != source['weights_sha256']:
        raise ValueError('Original failed fit changed')
    saved = np.load(wp); g = particle_geometry(root, '10028', 'inference_half0', radius=12, count=128, seed=609315)
    result = dict(complete=False, planned_rank=8192, planned_iterations=300, tau=2.,
        scope='Numerical only. Residual diagonal at rank2048 is ~9860 versus ridge .25; plan reuse alone improves matvec time only ~8%. Probe a substantially larger low-rank inverse on the original first case. No new pixels, maps, coverage or parameter selection by scientific outcomes.',
        source_snapshot=source_snapshot(root, Path(__file__), ['src/fourier_splats/uq_cached_quadrature.py']),
        initialization_sha256=source['weights_sha256'])
    def save():
        out.write_text(json.dumps(result, indent=2)+'\n')
    save(); tick = time.perf_counter()
    gram = CachedQuadratureGram(QuadratureObservationGram(g['k'], g['ctf'], float(saved['noise_std']), order=80, preconditioner_rank=8192))
    result.update(setup_seconds=time.perf_counter()-tick, preconditioner_diagnostics=gram.preconditioner_diagnostics)
    save(); print('setup', result['setup_seconds'], result['preconditioner_diagnostics'], flush=True)
    a, _ = gram.target(source['centers'], [1], source['width_fraction_field'])
    system = LinearOperator(gram.shape, matvec=lambda x: gram.matvec(x)+.25*x, dtype=float)
    pre = gram.preconditioner(.25); count = [0]; trace = []; tick = time.perf_counter()
    def callback(w):
        count[0] += 1
        if count[0] % 10 == 0:
            residual = float(np.linalg.norm(system@w-a)/np.linalg.norm(a))
            trace.append([count[0], residual]); print(count[0], residual, flush=True)
    w, info = cg(system, a, x0=saved['weights'], M=pre, rtol=1e-10, atol=0., maxiter=300, callback=callback)
    residual = system@w-a
    np.savez_compressed(out.with_suffix('.npz'), weights=w, indices=g['indices'], noise_std=float(saved['noise_std']))
    result.update(complete=True, fit_seconds=time.perf_counter()-tick, cg_info=int(info),
        iterations=count[0], relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(a)), trace=trace,
        variance_gap_diagnostic=float(4*residual@(pre@residual)),
        arrays_sha256=hashlib.sha256(out.with_suffix('.npz').read_bytes()).hexdigest())
    save(); print('done', result['relative_residual'], result['variance_gap_diagnostic'], flush=True)


if __name__ == '__main__':
    main()
