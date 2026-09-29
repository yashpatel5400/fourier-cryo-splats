#!/usr/bin/env python3
"""Search feasible poses with a surrogate; validate bias in the full L2 cube.

Quadrature and float32 acceleration are used only to propose adversarial poses.
The final bias eliminates the continuous density ball with independent analytic
sinc integrals and evaluates the fixed constant-cell pilot exactly. Results are
feasible lower bounds, never certified global nonlinear maxima.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import torch
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous import cell_forward, continuous_residual_norm
from fourier_splats.uq_continuous_pose import cube_quadrature, pose_cell_forward, perturbed_geometry
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_nonlinear_stress import TorchPoseAdjoint
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'


def run(args, dataset, snapshot):
    source_path = BASE/'continuous-quadrature-optimized'/f'{dataset}.json'
    source = json.loads(source_path.read_text()); cfg = source['config']
    g = particle_geometry(ROOT, dataset, 'inference_half0', radius=cfg['frequency_radius'],
                          count=cfg['particles'], seed=cfg['seed'])
    checkpoint = np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
    op, pilot, _, noise = model(g, checkpoint, 24); pilot = op.expand(pilot)
    nominal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
    xyz, quadrature = cube_quadrature(args.order)
    indices = np.minimum(((xyz+.5)*24).astype(int), 23)
    pilot_at_nodes = pilot.reshape(24, 24, 24)[indices[:, 2], indices[:, 1], indices[:, 0]]*24**1.5
    fake_op = SimpleNamespace(op=SimpleNamespace(k=g['k'], box=1., transfer=g['ctf']/noise,
                                               n=len(g['k']), q=g['k'].shape[1]), xyz=xyz,
                              shape=(2*np.prod(g['k'].shape[:2]), len(xyz)))
    device = args.device; dtype = torch.float32 if device == 'mps' else torch.float64
    tensor = lambda x: torch.as_tensor(np.asarray(x).copy(), dtype=dtype, device=device)
    qw = tensor(quadrature); pilotq = tensor(pilot_at_nodes)
    out = BASE/args.output; out.mkdir(exist_ok=True)
    result = {'stage': 'development continuous-density feasible nonlinear adversaries',
              'dataset': dataset, 'config': vars(args), 'source_snapshot': snapshot,
              'source_fit_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'search_scope': 'Quadrature surrogate proposes poses; exact continuous density ball at each final pose.',
              'not_claimed': 'Global nonlinear maximization, validated arithmetic, or empirical pose calibration.',
              'records': []}
    rng = np.random.default_rng(args.seed+int(dataset))
    for name in args.targets.split(','):
        centers = [[0, 0, 0]] if name == 'center' else [[0, 0, .08], [0, 0, -.08]]
        signs = [1] if name == 'center' else [1, -1]
        ell = sum(sign*np.exp(-np.sum((xyz-center)**2, axis=1)/(2*args.width**2))
                  /(2*np.pi*args.width**2)**1.5 for center, sign in zip(np.asarray(centers), signs))
        ellq = tensor(ell)
        saved = np.load(BASE/'continuous-quadrature-optimized'/f'{dataset}-{name}-{args.width}-weights.npz')
        w = saved['weights']; np.testing.assert_array_equal(saved['indices'], g['indices'])
        for degrees in map(float, args.angles.split(',')):
            started = time.perf_counter(); angle = np.deg2rad(degrees); shift = .01/24
            adjoint = TorchPoseAdjoint(fake_op, g['q'], w, angle, shift, device=device, dtype=dtype)
            zero = torch.zeros((len(g['k']), 5), device=device, dtype=dtype)
            with torch.no_grad(): nominalq = adjoint(zero)
            candidates = []
            for sign in [-1., 1.]:
                for start in range(args.starts):
                    initial = np.zeros((len(g['k']), 5))
                    if start == 1:
                        initial[:, 0] = sign
                    elif start > 1:
                        initial = rng.normal(size=initial.shape)
                        initial /= np.linalg.norm(initial, axis=1)[:, None]
                    u = torch.nn.Parameter(tensor(initial)); optimizer = torch.optim.Adam([u], lr=args.lr)
                    best_value = -np.inf; best_pose = None; trace = []
                    for iteration in range(args.steps+1):
                        optimizer.zero_grad(set_to_none=True)
                        field = adjoint(u)
                        value = sign*torch.sum(qw*pilotq*(field-nominalq))+2*torch.sqrt(torch.sum(qw*(field-ellq)**2))
                        scalar = float(value.detach().cpu())
                        if scalar > best_value:
                            best_value = scalar; best_pose = u.detach().cpu().numpy().astype(float).copy()
                        if iteration % 25 == 0: trace.append({'iteration': iteration, 'surrogate_bias': scalar})
                        if iteration < args.steps:
                            (-value).backward(); optimizer.step()
                            with torch.no_grad():
                                u.div_(torch.clamp(torch.linalg.vector_norm(u, dim=1, keepdim=True), min=1.))
                    # Float32 projection can overshoot one slightly. Reproject in
                    # float64 before the independent physical forward evaluation.
                    best_pose /= np.maximum(1., np.linalg.norm(best_pose, axis=1, keepdims=True))
                    rotated, shifts = perturbed_geometry(g['k'], g['q'], best_pose, angle, shift)
                    exact = continuous_residual_norm(rotated, g['ctf']*np.exp(2j*np.pi*shifts),
                                                     w, noise, centers, signs, args.width)
                    expected = pose_cell_forward(g['k'], g['q'], g['ctf'], pilot, 24, noise, best_pose, angle, shift)
                    pilot_bias = float(w@(expected-nominal))
                    # Unpadded value is the numerical feasible bias estimate;
                    # the padding is an upper error allowance, not a lower proof.
                    norm_h = np.sqrt(max(0., exact['residual_norm2_unpadded']))
                    signed_bias = sign*pilot_bias+2*norm_h
                    feasible = abs(pilot_bias)+2*norm_h
                    if start == 0:
                        # Independent SciPy rotation + direct Fourier sum checks
                        # the accelerated proposal model at a small node subset.
                        ids = np.arange(0, len(xyz), max(1, len(xyz)//31))
                        k_scipy = g['k']@Rotation.from_rotvec(angle*best_pose[:, :3]).as_matrix()
                        weights = w.reshape(len(g['k']), -1); count = g['k'].shape[1]
                        coeff = (weights[:, :count]+1j*weights[:, count:])*g['ctf']/noise
                        direct = np.zeros(len(ids))
                        for i in range(len(g['k'])):
                            phase = 2j*np.pi*(k_scipy[i]@xyz[ids].T+shifts[i, :, None])
                            direct += (coeff[i]@np.exp(phase)).real
                        predicted = adjoint(tensor(best_pose)).detach().cpu().numpy()[ids]
                        error = float(np.linalg.norm(direct-predicted)/max(np.linalg.norm(direct), 1e-12))
                        if error > (2e-4 if device == 'mps' else 1e-10):
                            raise AssertionError('Accelerated proposal disagrees with independent Fourier field')
                    else: error = None
                    candidates.append({'sign': sign, 'start': start, 'surrogate_maximum': best_value,
                                       'exact_signed_bias_numerical': float(signed_bias),
                                       'exact_absolute_bias_numerical': float(feasible),
                                       'pilot_bias': pilot_bias, 'analytic_residual': exact,
                                       'maximum_pose_norm': float(np.linalg.norm(best_pose, axis=1).max()),
                                       'independent_field_relative_error': error, 'trace': trace})
                    np.savez(out/f'{dataset}-{name}-{degrees}-sign{sign}-start{start}.npz', poses=best_pose)
                    print(dataset, name, degrees, sign, start, 'surrogate', best_value, 'exact', feasible, flush=True)
            result['records'].append({'target': name, 'width': args.width, 'angle_degrees': degrees,
                                      'candidates': candidates,
                                      'largest_feasible_bias_numerical': max(c['exact_absolute_bias_numerical'] for c in candidates),
                                      'seconds': time.perf_counter()-started})
            (out/f'{dataset}.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--datasets', default='10028,10049,10076')
    p.add_argument('--targets', default='center,contrast'); p.add_argument('--width', type=float, default=.07)
    p.add_argument('--angles', default='.5,1,2'); p.add_argument('--order', type=int, default=12)
    p.add_argument('--starts', type=int, default=3); p.add_argument('--steps', type=int, default=150)
    p.add_argument('--lr', type=float, default=.05); p.add_argument('--device', default='mps', choices=['cpu', 'mps'])
    p.add_argument('--seed', type=int, default=609501); p.add_argument('--output', default='continuous-pose-adversaries')
    args = p.parse_args()
    snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py'])
    for dataset in args.datasets.split(','): run(args, dataset, snapshot)
