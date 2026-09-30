#!/usr/bin/env python3
"""Post-fit W1/W4 sensitivity checks; never retune weights or replace a certificate."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import finufft
import numpy as np
from scipy.special import logsumexp
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_cubic_design import DifferentiableCubicPoseFieldOperator
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
SOURCE = BASE/'cubic-weight-probe/10049-pilot_region_1-1.json'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--wait-hours', type=float, default=0.)
    args = parser.parse_args(); deadline = time.monotonic()+3600*args.wait_hours
    if 'tmp/nufft-openmp' not in finufft.__file__ or 'torch' in sys.modules:
        raise RuntimeError('Isolated CPU-only FINUFFT required')
    while True:
        try:
            case = json.loads(SOURCE.read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            case = {}
        if case.get('error'):
            raise RuntimeError(f'Original fit failed; preserve it: {case["error"]}')
        if case.get('complete'):
            break
        if time.monotonic() >= deadline:
            raise TimeoutError('Awaited final weights are not yet available')
        time.sleep(15)
    path = BASE/'audit-regressions/cubic-design-final-numerics.json'
    if path.exists():
        raise RuntimeError('Preserve prior diagnostics')
    result = {'complete': False, 'source_case_sha256': sha(SOURCE),
        'scope': 'Two-tolerance empirical sensitivity, not an operator-error upper, new probability certificate or weight-selection rule.',
        'eps_original': 1e-12, 'eps_check': 1e-14, 'diagnostic_seed': 650129,
        'source_snapshot': source_snapshot(ROOT, Path(__file__), ['research/uncertainty/CUBIC-DESIGN-NUMERICAL-CHECK.md'])}
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    save(); start = time.perf_counter()
    try:
        fit = case['fit']; wp = SOURCE.with_suffix('.npz')
        if sha(wp) != case['arrays_sha256']:
            raise ValueError('Final weights changed')
        result['source_arrays_sha256'] = sha(wp)
        saved = np.load(wp); w = saved['weights']; noise = float(saved['noise_std'])
        old = json.loads((ROOT/case['source_fit']).read_text())
        g = particle_geometry(ROOT, '10049', 'inference_half0', radius=12, count=128, seed=old['config']['seed'])
        np.testing.assert_array_equal(g['indices'], saved['indices'])
        gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=80, eps=1e-12, preconditioner_rank=0)
        gram.nthreads = 2; original = gram.matvec(w); gram.eps = 1e-14; checked = gram.matvec(w)
        difference = checked-original
        quadratic = float(abs(w@difference)); action = float(np.linalg.norm(difference))
        guard = fit['density_squared_roundoff_guard']
        result['nominal_gram'] = {'absolute_quadratic_change': quadratic, 'gram_action_change_norm': action,
            'relative_action_change': action/max(float(np.linalg.norm(checked)), 1e-300),
            'original_roundoff_only_squared_guard': guard,
            'quadratic_change_over_guard': quadratic/max(guard, 1e-300),
            'original_density_bias': fit['density_bias'],
            'heuristic_density_bias_after_adding_observed_change': float(2*np.sqrt((fit['density_bias']/2)**2+quadratic)),
            'empirical_extra_dual_defect_at_original_factor': float(4*action/fit['density_bias']),
            'interpretation': 'Observed differences can be added as heuristic sensitivity margins; neither tolerance is validated ground truth.'}
        op = DifferentiableCubicPoseFieldOperator(g['k'], g['q'], g['ctf'], w, noise, np.deg2rad(1.), .5/g['field_A'],
            order=64, nthreads=2, eps=1e-12, block_particles=16)
        op.denominators = saved['column_denominators'].copy()
        rng = np.random.default_rng(650129)
        u = rng.normal(size=op.shape[1]); u /= np.linalg.norm(u)
        v = rng.normal(size=op.shape[0]); v /= np.linalg.norm(v)
        original_forward = op.matvec(u); original_adjoint = op.rmatvec(v)
        op.eps = 1e-14; checked_forward = op.matvec(u); checked_adjoint = op.rmatvec(v)
        result['cubic_field'] = {'forward_absolute_difference': float(np.linalg.norm(checked_forward-original_forward)),
            'forward_relative_difference': float(np.linalg.norm(checked_forward-original_forward)/np.linalg.norm(checked_forward)),
            'adjoint_absolute_difference': float(np.linalg.norm(checked_adjoint-original_adjoint)),
            'adjoint_relative_difference': float(np.linalg.norm(checked_adjoint-original_adjoint)/np.linalg.norm(checked_adjoint)),
            'original_adjoint_residual': float(abs(v@original_forward-u@original_adjoint)),
            'check_adjoint_residual': float(abs(v@checked_forward-u@checked_adjoint)),
            'scope': 'One fixed diagnostic pair; not a spectral/operator-error bound.'}
        selected = fit['optimization_history'][fit['selected_evaluation']-1]
        singular = np.asarray(selected['ritz_singular_values']); ritz = float(np.max(singular)**2)
        rayleigh = fit['spectral_upper_certificate']['rayleigh_lower']; upper = fit['spectral_upper_certificate']['eigenvalue_upper']
        result['selected_spectral_guide'] = {'selected_evaluation': fit['selected_evaluation'],
            'ritz_eigenvalue': ritz, 'fresh_rayleigh_lower': rayleigh, 'fresh_upper': upper,
            'rayleigh_over_ritz': rayleigh/ritz, 'upper_over_ritz': upper/ritz,
            'interpretation': 'Ritz values are design guides; the fresh Rayleigh value is another lower bound and the upper is the probabilistic audit.'}
        smooth = fit['spectral_smoothing']; value = smooth*logsumexp(singular/smooth)
        mixing = np.exp((singular-value)/smooth)
        result['gap_contributions'] = {'nominal_eta_gradient_norm': float(4*fit['density_majorant_eta']*np.linalg.norm(w)/fit['density_bias']),
            'spectral_entropy_term_before_outer_norm': float(2*case['lifted_radius']*(value-mixing@singular)),
            'maximum_smoothing_additive_bound': fit['smoothing_additive_objective_bound'],
            'roundoff_guard': guard,
            'scope': 'Components explaining why a valid lower need not be tight; not a decomposition of the observed primal-dual gap.'}
        result.update(complete=True, seconds=time.perf_counter()-start); save()
        print(json.dumps({k: v for k, v in result.items() if k != 'source_snapshot'}, indent=2), flush=True)
    except Exception as exc:
        result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
