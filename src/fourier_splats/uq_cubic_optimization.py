"""Bounded cubic pose-aware design with a separately drawn spectral audit.

Ritz modes guide an approximate optimizer. They are never an upper bound.
Quadrature pads hold in real arithmetic; floating-point guards are heuristic.
"""
import time
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.sparse.linalg import LinearOperator, eigsh, ArpackNoConvergence
from scipy.stats import norm
from fourier_splats.uq_pose_optimization import amplitude_gradient, gram_spectral_power
from fourier_splats.uq_continuous_pose import polynomial_kernel_error
from fourier_splats.uq_random_spectral import gaussian_power_upper
from fourier_splats.uq_intervals import bias_aware_half_width_stable


class CubicDesignObjective:
    """Convex norm surrogate, evaluated using supplied spatial Ritz directions.

    Exact leading singular directions give a spectral-norm subgradient when
    smoothing=0. Approximate directions still define valid norm supports for
    the dual bound, but their values do not certify the spectral norm above.
    """
    def __init__(self, gram, operator, remainder, moments, centers, signs, width,
                 density_radius, pilot_norm, lifted_radius, z, smoothing=.002):
        self.gram, self.op, self.remainder = gram, operator, remainder
        self.moments = np.asarray(moments)
        self.B, self.P, self.L, self.z = map(float, (density_radius, pilot_norm, lifted_radius, z))
        self.smoothing = float(smoothing)
        if min(self.B, self.L, self.z) <= 0 or self.P < 0 or self.smoothing < 0:
            raise ValueError('Positive radii/critical value and nonnegative pilot norm/smoothing required')
        self.a, self.target_norm2 = gram.target(centers, signs, width)
        self.coefficients = operator.coefficient_bound_matrix()
        self.kernel_error = polynomial_kernel_error(operator.k, operator.order, degree=6)
        # |w'(G_Q-G)w| <= E (sum |t| pair_amplitude)^2 <= eta ||w||^2.
        # Therefore G_Q+eta I majorizes the continuous Gram in real arithmetic.
        self.eta = float(gram.kernel_error*np.sum(gram.transfer**2))
        self.n, self.nq = operator.n, operator.nq

    def density(self, weights, final_guard=False):
        gw = self.gram.matvec(weights)
        linear = float(weights@self.a); quadratic = float(weights@gw)
        squared = float(self.target_norm2-2*linear+quadratic+self.eta*(weights@weights))
        guard = 256*np.finfo(float).eps*(abs(self.target_norm2)+2*abs(linear)+abs(quadratic))
        if squared < -max(guard, 1e-12*abs(self.target_norm2)):
            raise FloatingPointError('Substantially negative majorized density residual')
        if squared <= 0 and not final_guard:
            raise FloatingPointError('Nondifferentiable/nonpositive density residual in design')
        residual = float(np.sqrt(max(0., squared+(guard if final_guard else 0.))))
        return residual, gw, guard

    def mass_terms(self, weights):
        w = weights.reshape(self.n, 2*self.nq)
        amplitude = np.hypot(w[:, :self.nq], w[:, self.nq:])
        masses = np.einsum('naj,nj->na', self.coefficients, amplitude)
        gradient = amplitude_gradient(np.einsum('naj,na->nj', self.coefficients, masses), weights)
        return float(self.kernel_error*np.sum(masses**2)), self.kernel_error*gradient

    def evaluate(self, weights, modes):
        weights = np.asarray(weights, float)
        self.op.set_weights(weights)
        modes = np.asarray(modes, float)
        if modes.ndim != 2 or modes.shape[0] != self.op.shape[0] or not modes.shape[1]:
            raise ValueError('Nonempty spatial mode matrix required')
        lengths = np.linalg.norm(modes, axis=0)
        if not np.isfinite(lengths).all() or np.any(lengths == 0):
            raise ValueError('Finite nonzero spatial modes required')
        modes = modes/lengths
        adjoints = [self.op.rmatvec(v) for v in modes.T]
        singular = np.array([np.linalg.norm(u) for u in adjoints])
        if self.smoothing > 0:
            field = float(self.smoothing*logsumexp(singular/self.smoothing))
            mixing = np.exp((singular-field)/self.smoothing)
        else:
            field = float(singular.max()); mixing = np.eye(len(singular))[int(np.argmax(singular))]
        tiny = np.finfo(float).tiny
        field_gradient = sum(p*self.op.weight_gradient(u/max(s, tiny), v)
            for p, u, s, v in zip(mixing, adjoints, singular, modes.T))
        pad, pad_gradient = self.mass_terms(weights)
        pose_norm = float(np.sqrt(field**2+pad))
        pose_gradient = self.B*self.L*(field*field_gradient+pad_gradient)/max(pose_norm, tiny)
        pilot_vector = self.op.pair_moments(self.moments)
        pilot_norm = float(np.linalg.norm(pilot_vector))
        pilot_gradient = self.L*self.op.moment_weight_gradient(pilot_vector/max(pilot_norm, tiny), self.moments)
        rem, rem_gradient = self.remainder.value_gradient(weights)
        rem_gradient *= self.B+self.P
        residual, gw, _ = self.density(weights)
        sd = float(np.linalg.norm(weights))
        support = pose_gradient+pilot_gradient+rem_gradient
        gradient = self.z*weights/max(sd, tiny)+self.B*(gw+self.eta*weights-self.a)/residual+support
        terms = {'noise_term': self.z*sd, 'density_term': self.B*residual,
                 'pose_polynomial_term': self.B*self.L*pose_norm,
                 'pilot_polynomial_term': self.L*pilot_norm,
                 'quartic_term': (self.B+self.P)*rem}
        value = float(sum(terms.values()))
        row = dict(terms, surrogate_objective=value, gradient_norm=float(np.linalg.norm(gradient)),
                   ritz_singular_values=singular.tolist(), active_spectral_modes=int(np.sum(mixing > 1e-3)))
        # Even with log-sum-exp, these are supports of the unsmoothed norm sum:
        # the mixture is convex and field/sqrt(field^2+pad) <= 1. Entropy means
        # support@weights may be smaller than the smoothed value.
        return value, gradient, row, support

    def audit(self, weights, norm_support, spectral_upper):
        weights = np.asarray(weights, float); self.op.set_weights(weights)
        residual, gw, guard = self.density(weights, final_guard=True)
        pad, _ = self.mass_terms(weights)
        field = self.L*np.sqrt(spectral_upper+pad)
        pilot = self.L*np.linalg.norm(self.op.pair_moments(self.moments))
        rem, _ = self.remainder.value_gradient(weights)
        bias = self.B*(residual+field)+pilot+(self.B+self.P)*rem
        sd = float(np.linalg.norm(weights)); upper = self.z*sd+bias
        factor = self.B/max(residual, np.finfo(float).tiny)
        error = self.gram.quadrature_error(weights)
        # eta*w belongs to the majorant, not to the continuous A*h. Omitting it
        # here and padding G_Q's action preserves the actual density support.
        defect = float(np.linalg.norm(factor*(self.a-gw)-norm_support)+factor*error['gram_action_norm'])
        scale = min(1., self.z/max(defect, np.finfo(float).tiny))
        lower = max(0., float(scale*factor*(self.target_norm2-weights@self.a)))
        if lower > upper+1e-8*max(1., upper):
            raise FloatingPointError('Dual lower bound exceeds audited upper')
        return {'noise_sd': sd, 'bias': float(bias), 'density_bias': float(self.B*residual),
                'pose_field_upper': float(field), 'pose_polynomial_bias': float(self.B*field),
                'pilot_polynomial_bias': float(pilot), 'quartic_bias': float((self.B+self.P)*rem),
                'pose_integration_pad': pad, 'density_majorant_eta': self.eta,
                'density_squared_roundoff_guard': float(guard), 'target_norm': float(np.sqrt(self.target_norm2)),
                'sum_objective_upper': float(upper), 'dual_lower_bound': lower,
                'dual_defect_upper': defect, 'dual_rescaling': scale,
                'relative_sum_gap': float((upper-lower)/max(upper, np.finfo(float).tiny))}


def optimize_cubic_weights(objective, initial_weights, alpha=.05/12, numerical_delta=1e-6/12,
                           maxiter=30, max_evaluations=40, modes_count=3, ritz_tolerance=1e-3,
                           ritz_subspace=13, optimization_seed=650011, certificate_seed=650001,
                           power_probes=4, power_iterations=40, precondition=True,
                           callback=None, checkpoint_callback=None):
    if not 0 < numerical_delta < alpha < .5:
        raise ValueError('0 < numerical delta < alpha < .5 required')
    if abs(objective.z-norm.isf((alpha-numerical_delta)/2)) > 1e-12:
        raise ValueError('Objective critical value and final error budget disagree')
    if not 0 < modes_count < ritz_subspace <= objective.op.shape[0]:
        raise ValueError('Invalid Ritz subspace')
    start = time.perf_counter(); initial = np.asarray(initial_weights, float)
    if not np.isfinite(initial).all() or np.linalg.norm(initial) == 0:
        raise ValueError('Finite nonzero initial weights required')
    initial_residual, _, _ = objective.density(initial)
    ridge = float(objective.z*initial_residual/(objective.B*np.linalg.norm(initial)))
    gram = objective.gram
    def transform(v, power):
        return gram_spectral_power(gram, v, ridge, power) if precondition else v
    coordinates = transform(initial, .5)
    np.testing.assert_allclose(transform(coordinates, -.5), initial, rtol=1e-5, atol=1e-9)
    mode = np.random.default_rng(optimization_seed).normal(size=objective.op.shape[0]); mode /= np.linalg.norm(mode)
    history = []; best = {'value': np.inf}
    def evaluate(x):
        nonlocal mode
        weights = transform(x, -.5); objective.op.set_weights(weights)
        operator = LinearOperator((objective.op.shape[0],)*2, matvec=objective.op.spatial_gram, dtype=float)
        status = {'converged': True, 'requested_modes': modes_count, 'returned_modes': modes_count}
        evaluation_start = time.perf_counter()
        try:
            values, vectors = eigsh(operator, k=modes_count, which='LA', v0=mode,
                tol=ritz_tolerance, ncv=ritz_subspace, maxiter=100)
        except ArpackNoConvergence as exc:
            if exc.eigenvalues is None or not len(exc.eigenvalues):
                raise
            values, vectors = exc.eigenvalues, exc.eigenvectors
            status.update(converged=False, returned_modes=len(values), message=str(exc))
        vectors = vectors[:, np.argsort(values)[::-1]]
        vectors /= np.linalg.norm(vectors, axis=0); mode = vectors[:, 0]
        value, gradient, row, support = objective.evaluate(weights, vectors)
        row.update(evaluation=len(history)+1, eigensolver_status=status,
                   seconds=time.perf_counter()-evaluation_start, elapsed_seconds=time.perf_counter()-start)
        history.append(row)
        if value < best['value']:
            best.update(value=value, weights=weights.copy(), support=support.copy(), evaluation=len(history))
            if checkpoint_callback:
                checkpoint_callback(weights.copy(), row)
        if callback:
            callback(row)
        return value, transform(gradient, -.5)
    optimized = minimize(evaluate, coordinates, jac=True, method='L-BFGS-B',
        options={'maxiter': maxiter, 'maxfun': max_evaluations, 'maxls': 10, 'ftol': 1e-5, 'gtol': 1e-6, 'maxcor': 15})
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
        spectral_upper_certificate=spectral, preconditioner_ridge=ridge, coordinate_preconditioning=precondition,
        partial_mode_evaluations=sum(not r['eigensolver_status']['converged'] for r in history),
        maximum_evaluations_requested=max_evaluations, evaluation_budget_overshoot=max(0, len(history)-max_evaluations),
        evaluation_budget_scope='The library checks maxfun after a line search, with maxls=10; every evaluation retained.',
        alpha_total=alpha, alpha_noise=alpha-numerical_delta, alpha_numerical=numerical_delta,
        pose_quadrature_order=objective.op.order, density_quadrature_order=gram.order,
        spectral_smoothing=objective.smoothing, spectral_modes=modes_count, ritz_subspace=ritz_subspace,
        smoothing_additive_objective_bound=float(objective.B*objective.L*objective.smoothing*np.log(modes_count)),
        gap_scope='Unsmoothed order-matched triangle surrogate with density majorant; not the folded or joint-refined optimum.',
        numerical_scope='Analytic quadrature bounds in real arithmetic; heuristic floating-point guards, not interval arithmetic.')
    return audit
