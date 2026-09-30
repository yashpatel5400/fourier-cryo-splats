#!/usr/bin/env python3
"""Isolated CPU-only OpenMP timing; never imports the incompatible Torch runtime."""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
import finufft
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def main():
    path = BASE/'nufft-openmp-profile.json'
    if path.exists(): raise RuntimeError('Preserve previous timing results')
    binary = Path(finufft.__file__).parent/'libfinufft.dylib'
    if 'tmp/nufft-openmp' not in str(binary.resolve()):
        raise RuntimeError('Use the isolated CPU wheel through PYTHONPATH')
    if 'torch' in sys.modules: raise RuntimeError('Do not mix OpenMP runtimes')
    snapshot = source_snapshot(ROOT, Path(__file__))
    g = particle_geometry(ROOT, '10028', 'inference_half0', radius=5, count=128, seed=609315)
    saved = np.load(BASE/'continuous-quadrature-optimized/10028-center-0.07-weights.npz')
    np.testing.assert_array_equal(g['indices'], saved['indices'])
    op = PolynomialPoseFieldOperator(g['k'], g['q'], g['ctf'], saved['weights'],
                                     float(saved['noise_std']), np.deg2rad(1), .5/g['field_A'], order=32)
    op.establish_coefficient_scaling(); rng = np.random.default_rng(609881); records = []
    for iteration in range(4):
        v = rng.normal(size=op.shape[0]); v /= np.linalg.norm(v)
        outputs = {}; times = {}
        for threads in rng.permutation([1, 2, 4, 6]):
            op.nthreads = int(threads); start = time.perf_counter()
            outputs[int(threads)] = op.spatial_gram(v); times[int(threads)] = time.perf_counter()-start
        errors = {t: float(np.linalg.norm(y-outputs[1])/np.linalg.norm(outputs[1])) for t, y in outputs.items()}
        row = {'iteration': iteration, 'seconds': times, 'relative_discrepancy_from_one_thread': errors}
        records.append(row); print(row, flush=True)
    medians = {str(t): float(np.median([r['seconds'][t] for r in records])) for t in [1, 2, 4, 6]}
    result = {'complete': True, 'source_snapshot': snapshot, 'finufft_binary_path': str(binary.relative_to(ROOT)),
              'finufft_binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
              'torch_imported': 'torch' in sys.modules, 'records': records, 'median_seconds': medians,
              'speedup_relative_to_one_thread': {t: medians['1']/v for t, v in medians.items()}}
    path.write_text(json.dumps(result, indent=2)+'\n'); print(result['speedup_relative_to_one_thread'], flush=True)


if __name__ == '__main__': main()
