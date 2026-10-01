#!/usr/bin/env python3
"""Pure numerical speed/accuracy check of a fixed high-band Gram."""
from pathlib import Path
import json
import time
import numpy as np
from fourier_splats.uq_cached_quadrature import CachedQuadratureGram
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_provenance import source_snapshot


def main():
    root = Path(__file__).resolve().parents[1]
    out = root/'results/uncertainty/development/audit-regressions/cached-quadrature-runtime.json'
    if out.exists():
        raise RuntimeError('Preserve previous preflight')
    g = particle_geometry(root, '10028', 'inference_half0', radius=12, count=128, seed=609315)
    source = json.loads((root/'results/uncertainty/development/pilot-selected-fixed-v1/10028-pilot_region_1.json').read_text())
    base = QuadratureObservationGram(g['k'], g['ctf'], source['noise_std'], order=80, preconditioner_rank=0)
    rng = np.random.default_rng(927); w = rng.normal(size=base.shape[0])
    result = dict(complete=False, scope='No-pixels/no-coverage numerical transform comparison. FINUFFT plan reuse only.',
        source_snapshot=source_snapshot(root, Path(__file__), ['src/fourier_splats/uq_cached_quadrature.py', 'tests/test_cached_quadrature.py']), records=[])
    tick = time.perf_counter(); cached = CachedQuadratureGram(base)
    result['plan_setup_seconds'] = time.perf_counter()-tick
    expected = base.matvec(w)
    for name, gram in [('uncached', base), ('cached', cached)]:
        timings = []
        for _ in range(10):
            start = time.perf_counter(); actual = gram.matvec(w); timings.append(time.perf_counter()-start)
        result['records'].append(dict(implementation=name, seconds=timings, median_seconds=float(np.median(timings)),
            relative_difference=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected))))
    result['complete'] = True
    out.write_text(json.dumps(result, indent=2)+'\n')
    print([(r['implementation'], r['median_seconds'], r['relative_difference']) for r in result['records']], flush=True)


if __name__ == '__main__':
    main()
