"""Cubic design in fixed invertible coordinates, outside frozen study code.

This changes optimizer coordinates only. The forward and gradient maps differ
for nonsymmetric whitening. Original empirical fits keep their old solver.
No new empirical fit is authorized or executed by this module alone.
"""
from numbers import Integral
import time
import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator, eigsh, ArpackNoConvergence
from scipy.stats import norm
from .uq_random_spectral import gaussian_power_upper
from .uq_intervals import bias_aware_half_width_stable


def optimize_with_coordinates(objective, initial_weights, *, coordinates, certificate_seed,
                           optimization_seed, alpha=.05/12, numerical_delta=1e-6/12,
                           maxiter=30, max_evaluations=40, modes_count=3, ritz_tolerance=1e-3,
                           ritz_subspace=13, ritz_maxiter=100,
                           power_probes=4, power_iterations=40, maxls=10,
                           ftol=1e-5, gtol=1e-6,
                           callback=None, checkpoint_callback=None):
    for seed in [certificate_seed, optimization_seed]:
        if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
            raise ValueError('Explicit nonnegative integer random seeds required')
    if certificate_seed == optimization_seed:
        raise ValueError('Separate design and certificate streams required')
    if min(maxiter, max_evaluations, ritz_maxiter, maxls) < 1 or min(ftol, gtol) < 0:
        raise ValueError('Positive iteration budgets and nonnegative tolerances required')
    if not 0 < numerical_delta < alpha < .5:
        raise ValueError('0 < numerical delta < alpha < .5 required')
    if abs(objective.z-norm.isf((alpha-numerical_delta)/2)) > 1e-12:
        raise ValueError('Objective critical value and final error budget disagree')
    if not 0 < modes_count < ritz_subspace <= objective.op.shape[0]:
        raise ValueError('Invalid Ritz subspace')
    start = time.perf_counter(); initial = np.asarray(initial_weights, float)
    if not np.isfinite(initial).all() or np.linalg.norm(initial) == 0:
        raise ValueError('Finite nonzero initial weights required')
    gram = objective.gram
    initial_coordinates = coordinates.to_coordinates(initial)
    np.testing.assert_allclose(coordinates.to_weights(initial_coordinates), initial, rtol=1e-5, atol=1e-9)
    mode = np.random.default_rng(optimization_seed).normal(size=objective.op.shape[0]); mode /= np.linalg.norm(mode)
    history = []; best = {'value': np.inf}
    def evaluate(x):
        nonlocal mode
        weights = coordinates.to_weights(x); objective.op.set_weights(weights)
        operator = LinearOperator((objective.op.shape[0],)*2, matvec=objective.op.spatial_gram, dtype=float)
        status = {'converged': True, 'requested_modes': modes_count, 'returned_modes': modes_count}
        evaluation_start = time.perf_counter()
        try:
            values, vectors = eigsh(operator, k=modes_count, which='LA', v0=mode,
                tol=ritz_tolerance, ncv=ritz_subspace, maxiter=ritz_maxiter)
        except ArpackNoConvergence as exc:
            if exc.eigenvalues is None or not len(exc.eigenvalues):
                raise
            values, vectors = exc.eigenvalues, exc.eigenvectors
            status.update(converged=False, returned_modes=len(values), message=str(exc))
        vectors = vectors[:, np.argsort(values)[::-1]]
        vectors /= np.linalg.norm(vectors, axis=0); mode = vectors[:, 0]
        value, gradient, row, support = objective.evaluate(weights, vectors)
        coordinate_gradient = coordinates.gradient_to_coordinates(gradient)
        if not (np.isfinite(value) and np.isfinite(coordinate_gradient).all()):
            raise FloatingPointError('Nonfinite objective or coordinate gradient')
        row.update(coordinate_gradient_norm=float(np.linalg.norm(coordinate_gradient)), evaluation=len(history)+1, eigensolver_status=status,
                   seconds=time.perf_counter()-evaluation_start, elapsed_seconds=time.perf_counter()-start)
        history.append(row)
        if value < best['value']:
            best.update(value=value, weights=weights.copy(), support=support.copy(), evaluation=len(history))
            if checkpoint_callback:
                checkpoint_callback(weights.copy(), row)
        if callback:
            callback(row)
        return value, coordinate_gradient
    optimized = minimize(evaluate, initial_coordinates, jac=True, method='L-BFGS-B',
        options={'maxiter': maxiter, 'maxfun': max_evaluations, 'maxls': maxls, 'ftol': ftol, 'gtol': gtol, 'maxcor': 15})
    weights = best['weights']; objective.op.set_weights(weights)
    optimization_seconds = time.perf_counter()-start
    if callback:
        callback({'stage': 'fresh_spectral_audit', 'evaluations': len(history), 'selected_evaluation': best['evaluation']})
    spectral = gaussian_power_upper(objective.op.spatial_gram, objective.op.shape[0], numerical_delta,
        power_probes, power_iterations, certificate_seed, callback=callback)
    audit = objective.audit(weights, best['support'], spectral['eigenvalue_upper'])
    audit.update(weights=weights, half_width=bias_aware_half_width_stable(audit['noise_sd'], audit['bias'], alpha-numerical_delta),
        optimizer_success=bool(optimized.success), optimizer_message=str(optimized.message),
        optimizer_iterations=int(optimized.nit), optimization_history=history,
        selected_evaluation=best['evaluation'], selected_approximate_objective=best['value'],
        optimization_seconds=optimization_seconds, certification_seconds=time.perf_counter()-start-optimization_seconds,
        spectral_upper_certificate=spectral, coordinate_metric_diagnostics=coordinates.diagnostics,
        optimization_seed=int(optimization_seed), certificate_seed=int(certificate_seed),
        design_options=dict(maxiter=maxiter, maxfun=max_evaluations, maxls=maxls, ftol=ftol, gtol=gtol, ritz_maxiter=ritz_maxiter),
        partial_mode_evaluations=sum(not r['eigensolver_status']['converged'] for r in history),
        maximum_evaluations_requested=max_evaluations, evaluation_budget_overshoot=max(0, len(history)-max_evaluations),
        evaluation_budget_scope='The library checks maxfun after a line search; maxls is recorded and every evaluation retained.',
        alpha_total=alpha, alpha_noise=alpha-numerical_delta, alpha_numerical=numerical_delta,
        pose_quadrature_order=objective.op.order, density_quadrature_order=gram.order,
        spectral_smoothing=objective.smoothing, spectral_modes=modes_count, ritz_subspace=ritz_subspace,
        smoothing_additive_objective_bound=float(objective.B*objective.L*objective.smoothing*np.log(modes_count)),
        gap_scope='Unsmoothed order-matched triangle surrogate with density majorant; not the folded or joint-refined optimum.',
        numerical_scope='Analytic quadrature bounds in real arithmetic; heuristic floating-point guards, not interval arithmetic.')
    return audit
