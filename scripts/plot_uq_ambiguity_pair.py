#!/usr/bin/env python3
"""Visualize cell averages of a continuous testing witness, not a molecular map."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous import cell_adjoint, cell_target_coefficients
from fourier_splats.uq_continuous_pose import perturbed_geometry
from fourier_splats.uq_cell_pose_pair import refine_cell_density
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def main():
    p = argparse.ArgumentParser(); p.add_argument('--source', required=True); p.add_argument('--name', required=True); args = p.parse_args()
    source = BASE/args.source; record = json.loads(source.read_text()); array_path = source.with_suffix('.npz')
    if not record.get('complete') or record.get('error'):
        raise ValueError('Completed witness required')
    if hashlib.sha256(array_path.read_bytes()).hexdigest() != record['witness_array_sha256']:
        raise ValueError('Witness changed')
    output = ROOT/'paper/figures'/args.name
    if output.with_suffix('.pdf').exists():
        raise RuntimeError('Do not overwrite earlier figure')
    dataset = record['dataset']; target = record['target']; box = 48; noise = record['noise_std']
    g = particle_geometry(ROOT, dataset, 'inference_half0', radius=5, count=128, seed=record['source_geometry_config']['seed'])
    checkpoint = BASE/'representation'/dataset/'real_particles-spacing-2.0.npz'
    op, pilot, _, _ = model(g, np.load(checkpoint), 24); pilot = refine_cell_density(op.expand(pilot), 24, box)
    centers = [[0, 0, 0]] if target == 'center' else [[0, 0, .08], [0, 0, -.08]]
    signs = [1] if target == 'center' else [1, -1]
    ell = cell_target_coefficients(box, centers, signs, record['width_fraction_field'])
    arrays = np.load(array_path); n, nq = g['ctf'].shape; w = arrays['weights'].reshape(n, 2*nq)
    angle = np.deg2rad(record['rotation_radius_degrees']); shift = record['translation_radius_A']/g['field_A']; fields = []
    for j in range(2):
        k, phase = perturbed_geometry(g['k'], g['q'], arrays[f'pose{j}'], angle, shift)
        rw = (w[:, :nq]+1j*w[:, nq:])*np.exp(2j*np.pi*phase)
        a = cell_adjoint(k, g['ctf'], np.concatenate([rw.real, rw.imag], axis=1).ravel(), box, noise)
        density = float(arrays['common_density_scale'])*(pilot+(-1 if j == 0 else 1)*float(arrays['residual_amplitudes'][j])*(ell-a))
        fields.append(density.reshape((box,)*3)*box**1.5)
    # Midplane is the average of the two central cell slices; no clipping or mask.
    slices = [field[box//2-1:box//2+1].mean(axis=0) for field in fields]
    images = [*slices, slices[1]-slices[0]]; limit = max(float(np.max(abs(im))) for im in images)
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.0), constrained_layout=True)
    extent = np.array([-1, 1, -1, 1])*g['field_A']/2
    for axis, data, title in zip(axes, images, [r'Density $\rho_0$', r'Density $\rho_1$', r'Difference $\rho_1-\rho_0$']):
        im = axis.imshow(data, origin='lower', extent=extent, cmap='RdBu_r', vmin=-limit, vmax=limit, interpolation='nearest')
        axis.set_title(title); axis.set_xlabel(r'$x$ ($\AA$)'); axis.set_ylabel(r'$y$ ($\AA$)')
    fig.colorbar(im, ax=axes, shrink=.78, label='Signed normalized density')
    fig.savefig(output.with_suffix('.pdf')); fig.savefig(output.with_suffix('.png'), dpi=180); plt.close(fig)
    metadata = {'source': str(source.relative_to(ROOT)), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                'source_array_sha256': record['witness_array_sha256'], 'display_cell_box': box,
                'display': 'Exact cell-average projections of continuous witness; averaged two central slices; common symmetric color scale, no clipping or mask',
                'scope': 'Illustration of signed mathematical class; not an experimental map or biologically admissible pair'}
    output.with_suffix('.json').write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    main()
