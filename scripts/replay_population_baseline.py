#!/usr/bin/env python3
"""Replay unchanged population baseline functions; independently solve two-state fits."""
import contextlib
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import pickle
import platform
import subprocess
import time
from pathlib import Path

import jax
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
AUTHOR = ROOT / 'background/counting_particles_paper'
COMMIT = 'b9099e1b7fb3f03207f94d89b153839bfcd6a7c4'
OUT = ROOT / 'results/uncertainty/development/population-author-replay-v1'


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def numpy_from_jax_pickle(fun, args, state, aval):
    if fun is not np._core.multiarray._reconstruct:
        raise ValueError('Unexpected reconstruction callable')
    value = fun(*args)
    value.__setstate__(state)
    if value.dtype.hasobject:
        raise ValueError('Object arrays are not accepted')
    return value


class ArrayOnlyUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        allowed = {
            ('numpy', 'ndarray'): np.ndarray,
            ('numpy', 'dtype'): np.dtype,
            ('numpy.core.multiarray', '_reconstruct'): np._core.multiarray._reconstruct,
            ('numpy._core.multiarray', '_reconstruct'): np._core.multiarray._reconstruct,
            ('numpy.core.multiarray', 'scalar'): np._core.multiarray.scalar,
            ('numpy._core.multiarray', 'scalar'): np._core.multiarray.scalar,
            ('jax._src.array', '_reconstruct_array'): numpy_from_jax_pickle,
        }
        if (module, name) not in allowed:
            raise ValueError(f'Unexpected pickle reference: {module}.{name}')
        return allowed[module, name]


def load(path):
    with path.open('rb') as f:
        return ArrayOnlyUnpickler(f).load()


def solve_two(log_likelihood):
    ell = np.asarray(log_likelihood, dtype=np.float64)
    if ell.ndim != 2 or ell.shape[1] != 2 or not np.isfinite(ell).all():
        raise ValueError('Finite two-column likelihoods required')
    likelihood = np.exp(ell - ell.max(axis=1, keepdims=True))
    a, b = likelihood.T
    difference = a - b

    def derivative(t):
        with np.errstate(divide='ignore', invalid='raise'):
            return np.mean(difference / (t * a + (1 - t) * b))

    if np.all(difference == 0):
        optimum, identified = .5, False
    elif derivative(0) <= 0:
        optimum, identified = 0., True
    elif derivative(1) >= 0:
        optimum, identified = 1., True
    else:
        optimum = brentq(derivative, 0., 1., xtol=1e-14)
        identified = True
    return optimum, likelihood, identified


def main():
    if OUT.exists():
        raise ValueError('Preserve every attempted replay')
    actual = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=AUTHOR, text=True).strip()
    if actual != COMMIT or subprocess.check_output(['git', 'status', '--porcelain'], cwd=AUTHOR).strip():
        raise ValueError('Author source must be the unchanged pinned checkout')
    for rel in ['scripts/replay_population_baseline.py', 'research/uncertainty/POPULATION-BASELINE-REPLAY-PROTOCOL.md']:
        if subprocess.check_output(['git', 'show', 'HEAD:' + rel], cwd=ROOT) != (ROOT / rel).read_bytes():
            raise ValueError('Commit protocol and runner before execution')
    jax.config.update('jax_enable_x64', False)
    spec = importlib.util.spec_from_file_location('counting_author_utils', AUTHOR / 'utils.py')
    author = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(author)
    OUT.mkdir(parents=True)
    started = time.perf_counter()
    result = dict(complete=False, author_commit=actual, git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  python=platform.python_version(), versions={n: importlib.metadata.version(n) for n in ['numpy', 'scipy', 'jax', 'jaxlib', 'cvxpy']},
                  jax_devices=[str(x) for x in jax.devices()], jax_enable_x64=False, inputs={}, cases=[])

    def save():
        (OUT / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')

    files = ['utils.py', 'LICENSE', 'pyproject.toml', 'compute_weights_spike_synthetic.py', 'compute_weights_spike_experimental.py']
    files += [str(p.relative_to(AUTHOR)) for p in (AUTHOR / 'data').rglob('*') if p.is_file()]
    result['inputs'] = {p: dict(sha256=sha(AUTHOR / p), bytes=(AUTHOR / p).stat().st_size) for p in sorted(files)}
    save()
    paths = [('synthetic', AUTHOR / f'data/spike_synthetic/likelihoods_assignments_dataset{i}.pkl') for i in range(10)]
    paths += [('experimental', p) for p in sorted((AUTHOR / 'data/spike_experimental').glob('log_likelihoods_*.npy'))]
    for kind, path in paths:
        begin = time.perf_counter()
        if kind == 'synthetic':
            data = load(path)
            ell, hard, error = data['log_likelihoods'], data['hard_assignments'], float(data['error_predicted'])
            truth = np.bincount(data['true_assignments'], minlength=2) / len(data['true_assignments'])
            tolerance = 1e-8
        else:
            ell = np.load(path, allow_pickle=False)
            hard = np.argmax(ell, axis=1)
            error = float(np.load(path.with_name(path.name.replace('log_likelihoods_', 'error_predicted_')), allow_pickle=False).item())
            truth, tolerance = None, 1e-3
        if not np.array_equal(hard, np.argmax(ell, axis=1)):
            raise ValueError('Saved hard assignments disagree with likelihood argmax')
        text = io.StringIO()
        with contextlib.redirect_stdout(text):
            weights = np.asarray(author.multiplicative_gradient(ell, tol=tolerance), dtype=float)
            soft = np.asarray(author.multiplicative_gradient(ell, max_iterations=1, tol=-1), dtype=float)
            deconvolved = np.asarray(author.deconvolve_assignments(hard, error), dtype=float)
        stdout = text.getvalue()
        (OUT / (kind + '-' + path.stem + '.stdout.txt')).write_text(stdout)
        opt, likelihood, identified = solve_two(ell)
        observed = float(np.mean(hard == 0))
        analytic_deconv = .5 if error == .5 else float(np.clip((observed - error) / (1 - 2 * error), 0, 1))
        denom = likelihood @ weights
        gap = float(np.max(np.mean(likelihood / denom[:, None], axis=0)) - 1)
        objective_loss = float(np.mean(np.log(opt * likelihood[:, 0] + (1 - opt) * likelihood[:, 1]) - np.log(denom)))
        row = dict(kind=kind, input=str(path.relative_to(AUTHOR)), particles=len(ell), tolerance=tolerance,
                   original_exit_message=stdout.count('exiting!') > 0,
                   author_weights=weights.tolist(), independent_mle=[opt, 1-opt], identified=identified,
                   maximum_weight_difference=float(np.max(abs(weights - [opt, 1-opt]))),
                   float64_simplex_gradient_gap=gap, mean_objective_loss=objective_loss,
                   hard_weights=[observed, 1-observed], soft_weights=soft.tolist(), deconvolved=deconvolved.tolist(),
                   analytic_deconvolved=[analytic_deconv, 1-analytic_deconv],
                   maximum_deconvolution_difference=float(np.max(abs(deconvolved - [analytic_deconv, 1-analytic_deconv]))),
                   true_empirical_fractions=None if truth is None else truth.tolist(),
                   predicted_misclassification=error, seconds=time.perf_counter()-begin)
        result['cases'].append(row)
        save()
        print(json.dumps(row), flush=True)
    if len(result['cases']) != 23:
        raise ValueError('Unexpected published-case inventory')
    if subprocess.check_output(['git', 'status', '--porcelain'], cwd=AUTHOR).strip():
        raise ValueError('Author checkout changed')
    result.update(complete=True, seconds=time.perf_counter()-started)
    save()


if __name__ == '__main__':
    main()
