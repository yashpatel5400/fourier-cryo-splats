#!/usr/bin/env python3
"""Bounded native RELION VDAM execution check; no structural accuracy claim."""
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time

import mrcfile
import numpy as np
import pandas as pd
import starfile
from fourier_splats.physics import ctf
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/uncertainty/development/relion-runtime-probe'
DATA = ROOT/'data/uncertainty/relion-runtime-probe'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if OUT.exists() or DATA.exists():
        raise RuntimeError('Preserve previous input and output records')
    OUT.mkdir(parents=True); DATA.mkdir(parents=True)
    exe = ROOT/'tmp/relion-env/bin/relion_refine'
    source = ROOT/'data/10028'
    manifest = json.loads((source/'manifest.json').read_text())
    meta = np.load(source/'metadata.npz')
    rows = list(csv.DictReader((ROOT/'research/uncertainty/splits/10028.csv').open()))
    pool = [int(r['output_index']) for r in rows if r['split'] == 'pilot']
    ids = np.random.default_rng(670001).choice(pool, 256, replace=False)
    images = np.asarray(np.load(source/'images.npy', mmap_mode='r')[ids], dtype=float)
    d = images.shape[-1]; pixel = manifest['pixel_size_A']; field = d*pixel
    yy, xx = np.meshgrid(np.arange(d)-d//2, np.arange(d)-d//2, indexing='ij')
    background = xx**2+yy**2 > (.43*d)**2
    mean = images[:, background].mean(axis=1)
    scale = images[:, background].std(axis=1)
    if not np.isfinite(images).all() or np.min(scale) <= 0:
        raise ValueError('Invalid image normalization inputs')
    normalized = (-manifest['data_sign']*(images-mean[:, None, None])
                  /scale[:, None, None]).astype(np.float32)
    with mrcfile.new(DATA/'particles.mrcs') as m:
        m.set_data(normalized); m.voxel_size = pixel
    params = meta['ctf'][ids].astype(float)
    optics_constants = params[:, [5, 6, 7]]
    if not np.all(optics_constants == optics_constants[0]):
        raise ValueError('Probe requires one homogeneous optics group')
    optics = pd.DataFrame({'rlnOpticsGroupName': ['opticsGroup1'], 'rlnOpticsGroup': [1],
        'rlnVoltage': [params[0, 5]], 'rlnSphericalAberration': [params[0, 6]],
        'rlnAmplitudeContrast': [params[0, 7]], 'rlnImagePixelSize': [pixel],
        'rlnImageSize': [d], 'rlnImageDimensionality': [2]})
    particles = pd.DataFrame({'rlnImageName': [f'{i+1:06d}@particles.mrcs' for i in range(len(ids))],
        'rlnOpticsGroup': np.ones(len(ids), dtype=int), 'rlnDefocusU': params[:, 2],
        'rlnDefocusV': params[:, 3], 'rlnDefocusAngle': params[:, 4],
        'rlnPhaseShift': params[:, 8]})
    starfile.write({'optics': optics, 'particles': particles}, DATA/'particles.star')
    q = np.column_stack([xx.ravel(), yy.ravel()])/field
    r2 = np.sum(q*q, axis=1)[None, :]; az = np.arctan2(q[:, 1], q[:, 0])[None, :]
    df = .5*(params[:, 2, None]+params[:, 3, None]+(params[:, 2, None]-params[:, 3, None])
              *np.cos(2*(az-np.deg2rad(params[:, 4, None]))))
    voltage = params[:, 5, None]*1000
    wavelength = 12.2643247/np.sqrt(voltage*(1+voltage*.978466e-6))
    phase = -np.pi*wavelength*df*r2+.5*np.pi*params[:, 6, None]*1e7*wavelength**3*r2**2
    gamma = phase-np.deg2rad(params[:, 8, None])-np.arcsin(params[:, 7, None])
    convention_error = float(np.max(np.abs(-np.sin(gamma)+ctf(q, params))))
    if convention_error > .001:
        raise ValueError(f'Unexpected CTF convention discrepancy: {convention_error}')
    command = [str(exe), '--i', 'particles.star', '--o', str(OUT/'run'),
        '--grad', '--denovo_3dref', '--iter', '5', '--K', '1', '--sym', 'C1',
        '--ctf', '--flatten_solvent', '--zero_mask', '--dont_combine_weights_via_disc',
        '--preread_images', '--pool', '30', '--pad', '1', '--particle_diameter', str(.8*field),
        '--oversampling', '1', '--healpix_order', '1', '--offset_range', '6', '--offset_step', '2',
        '--auto_sampling', '--tau2_fudge', '4', '--j', '2', '--random_seed', '670001',
        '--grad_ini_subset', '128', '--grad_fin_subset', '256', '--grad_ini_resol', '60',
        '--grad_fin_resol', '25', '--grad_write_iter', '1']
    record = {'complete': False, 'stage': 'Native CPU runtime probe only; not a scientific reconstruction baseline',
        'dataset': '10028', 'particles': len(ids), 'indices': ids.tolist(), 'seed': 670001,
        'timeout_seconds': 600, 'command': command, 'cwd': str(DATA.relative_to(ROOT)),
        'version': subprocess.check_output([str(exe), '--version'], text=True, stderr=subprocess.STDOUT),
        'executable_sha256': sha(exe), 'input_hashes': {str(p.relative_to(ROOT)): sha(p) for p in
             [DATA/'particles.star', DATA/'particles.mrcs', source/'metadata.npz', source/'manifest.json']},
        'normalization': {'raw_pixel_multiplier': -manifest['data_sign'], 'background_radius_fraction': .43,
            'subtract_background_mean': True, 'scale_range': [float(scale.min()), float(scale.max())]},
        'ctf_sign_comparison_max_absolute_error': convention_error,
        'ctf_comparison_scope': 'Analytic RELION 5.0.1 formula versus project formula, retaining different wavelength constants; not binary operator conformance',
        'supplied_poses_or_maps': False,
        'source_snapshot': source_snapshot(ROOT, Path(__file__), ['research/uncertainty/RELION-RUNTIME-PROTOCOL.md'])}
    path = OUT/'profile.json'; path.write_text(json.dumps(record, indent=2)+'\n')
    log = ROOT/'logs/uncertainty/relion-abinit-profile.log'; started = time.perf_counter()
    try:
        with log.open('w') as f:
            result = subprocess.run(command, cwd=DATA, stdout=f, stderr=subprocess.STDOUT,
                timeout=600, env={**os.environ, 'OMP_NUM_THREADS': '2', 'OPENBLAS_NUM_THREADS': '1'})
        record.update(returncode=result.returncode, status='finished' if result.returncode == 0 else 'failed')
    except subprocess.TimeoutExpired:
        record.update(returncode=None, status='stopped at predeclared wall-clock ceiling')
    record.update(complete=True, seconds=time.perf_counter()-started,
        child_peak_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        log_sha256=sha(log), output_files={p.name: {'sha256': sha(p), 'bytes': p.stat().st_size}
            for p in OUT.iterdir() if p.is_file() and p != path})
    path.write_text(json.dumps(record, indent=2)+'\n')
    print({k: record[k] for k in ['complete','status','seconds','child_peak_rss_bytes']}, flush=True)


if __name__ == '__main__':
    main()
