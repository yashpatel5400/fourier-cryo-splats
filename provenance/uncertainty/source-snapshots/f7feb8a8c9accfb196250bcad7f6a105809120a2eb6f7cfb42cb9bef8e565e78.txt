#!/usr/bin/env python3
"""Bounded CPU timing probe of stock homogeneous cryoDRGN-AI pose search.

Thirty-two old development particles and two epochs cannot validate a recovered
structure. This records feasibility and timing only, with a 15-minute ceiling.
"""
import csv
import hashlib
import json
import os
import pickle
import subprocess
import time
from pathlib import Path
import cryodrgn
import numpy as np
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/uncertainty/development/abinit-runtime-probe'
OUT.mkdir(parents=True, exist_ok=True)
rows = list(csv.DictReader((ROOT/'research/uncertainty/splits/10028.csv').open()))
pool = np.array([int(r['output_index']) for r in rows if r['split'] == 'pilot'])
indices = np.random.default_rng(609491).choice(pool, 32, replace=False)
ind = ROOT/'data/uncertainty/abinit-runtime-probe-indices.pkl'
with ind.open('wb') as file:
    pickle.dump(indices, file)
source = ROOT/'data/10028/cryodrgn'
command = [str(ROOT/'.venv/bin/cryodrgn'), 'abinit', str(source/'particles.mrcs'),
           '--ctf', str(source/'ctf.pkl'), '--ind', str(ind), '-o', str(OUT/'run'),
           '--zdim', '0', '--num-epochs', '2', '--epochs-pose-search', '1',
           '--n-imgs-pretrain', '32', '--batch-size-hps', '2',
           '--batch-size-sgd', '8', '--num-workers', '0', '--max-threads', '4',
           '--no-amp', '--no-analysis', '--seed', '609491', '--log-interval', '8']
manifest = json.loads((ROOT/'data/10028/manifest.json').read_text())
if manifest['data_sign'] == 1:
    command.append('--uninvert-data')
result = {'stage': 'timing only, not converged reconstruction or scientific baseline',
          'dataset': '10028', 'particle_count': 32, 'timeout_seconds': 900,
          'stock_default_neural_architecture': '3 hidden layers, width 256, Hartley field',
          'stock_device_policy': 'CUDA if available, otherwise CPU; no MPS path',
          'input': 'Old development pilot particles; no supplied poses',
          'cryodrgn_version': cryodrgn.__version__, 'command': command,
          'indices_sha256': hashlib.sha256(ind.read_bytes()).hexdigest(),
          'source_snapshot': source_snapshot(ROOT, Path(__file__)), 'complete': False}
(OUT/'profile.json').write_text(json.dumps(result, indent=2)+'\n')
started = time.perf_counter()
try:
    with (ROOT/'logs/uncertainty/cryodrgn-abinit-profile.log').open('w') as log:
        run = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=900,
                             env={**os.environ, 'OMP_NUM_THREADS': '4', 'OPENBLAS_NUM_THREADS': '4'})
    result['returncode'] = run.returncode
    result['status'] = 'finished' if run.returncode == 0 else 'failed'
except subprocess.TimeoutExpired:
    result['returncode'] = None
    result['status'] = 'stopped at predeclared runtime ceiling'
result['seconds'] = time.perf_counter()-started
result['complete'] = True
(OUT/'profile.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2), flush=True)
