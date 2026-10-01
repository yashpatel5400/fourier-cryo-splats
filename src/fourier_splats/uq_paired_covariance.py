"""Matrix bilinear Gaussian ingredients; no continuous-pose certificate.

The inputs are real coordinates after independently fixed whitening and an
orthonormal projection. Each exposure covariance must be bounded by I in
Loewner order, and the two exposure errors must be independent.
"""
import numpy as np
from scipy.optimize import linprog, minimize_scalar


def _symmetric(matrix):
    a = np.asarray(matrix, float)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or not np.isfinite(a).all():
        raise ValueError('Finite square matrix required')
    if not np.allclose(a, a.T, rtol=1e-12, atol=1e-12):
        raise ValueError('Symmetric matrix required')
    return (a+a.T)/2


def matrix_gaussian_terms(weights):
    """Return the log normalizer and common-mean quadratic coefficient."""
    t = _symmetric(weights)
    eigenvalues, vectors = np.linalg.eigh(t)
    if np.max(np.abs(eigenvalues)) >= 1:
        raise ValueError('All eigenvalues must lie strictly between -1 and 1')
    normalizer = .5*np.log1p(-eigenvalues**2).sum()
    u = (vectors*(eigenvalues/(1-eigenvalues)))@vectors.T
    return float(normalizer), u


def matrix_log_factor(first, second, weights, null_exponent_upper):
    t = _symmetric(weights)
    x, z = np.asarray(first, float), np.asarray(second, float)
    if x.shape != (len(t),) or z.shape != x.shape or not np.isfinite([x,z]).all():
        raise ValueError('Two finite real observation vectors required')
    if not np.isfinite(null_exponent_upper) or null_exponent_upper < 0:
        raise ValueError('Externally justified nonnegative null-exponent bound required')
    normalizer, _ = matrix_gaussian_terms(t)
    return float(x@t@z+normalizer-null_exponent_upper)


def svec(matrix):
    a = np.asarray(matrix, float)
    i, j = np.triu_indices(a.shape[-1])
    return a[..., i, j]*np.where(i == j, 1., np.sqrt(2.))


def smat(vector, dimension):
    i, j = np.triu_indices(dimension)
    a = np.zeros((dimension, dimension))
    a[i, j] = np.asarray(vector)/np.where(i == j, 1., np.sqrt(2.))
    a[j, i] = a[i, j]
    return a


def covariance_cone_diagnostic(means, signal_second_moment, time_limit=60.):
    """Finite-view cone approximation and a numerically repaired separator.

    The LP minimizes a coordinate-scaled infinity residual. Its dual gives a
    separating quadratic form. Neither optimization nor floating-point repair
    establishes validity between the supplied views.
    """
    m = np.asarray(means, float); s = _symmetric(signal_second_moment)
    if m.ndim != 2 or m.shape[1] != len(s) or not np.isfinite(m).all():
        raise ValueError('Finite matching view means required')
    if np.linalg.eigvalsh(s).min() < -1e-10 or np.trace(s) <= 0:
        raise ValueError('Nonzero positive-semidefinite signal moment required')
    norms = np.sum(m*m, axis=1)
    nonzero = norms > 0
    if not nonzero.any():
        raise ValueError('At least one nonzero null view required')
    i, j = np.triu_indices(len(s))
    selected_means = m[nonzero]
    features = selected_means[:, i]*selected_means[:, j]*np.where(i == j, 1., np.sqrt(2.))
    features /= norms[nonzero, None]
    target = svec(s)
    scale = np.maximum(np.abs(target), 1e-3*np.max(np.abs(target)))
    a = features.T/scale[:, None]; b = target/scale
    constraints = np.vstack([np.column_stack([a, -np.ones(len(b))]),
                              np.column_stack([-a, -np.ones(len(b))])])
    objective = np.zeros(a.shape[1]+1); objective[-1] = 1
    result = linprog(objective, A_ub=constraints, b_ub=np.concatenate([b,-b]),
        bounds=(0, None), method='highs', options=dict(time_limit=float(time_limit),
            primal_feasibility_tolerance=1e-9, dual_feasibility_tolerance=1e-9))
    if not result.success:
        raise ArithmeticError(f'Cone LP status {result.status}: {result.message}')
    coefficients = np.zeros(len(m))
    coefficients[nonzero] = np.maximum(result.x[:-1], 0)/norms[nonzero]
    approximation = np.einsum('n,ni,nj->ij', coefficients, m, m)
    residual = np.max(np.abs(svec(approximation)-target)/scale)
    marginal = result.ineqlin.marginals
    dual_vector = (marginal[:len(b)]-marginal[len(b):])/scale
    raw_direction = smat(dual_vector, len(s))
    direction_norm = np.linalg.norm(raw_direction, 2)
    direction = raw_direction/max(direction_norm, 1e-30)
    violations = np.einsum('ni,ij,nj->n', m, direction, m)[nonzero]/norms[nonzero]
    correction = max(0., float(violations.max()))+1e-12
    direction -= correction*np.eye(len(s))
    repaired = np.einsum('ni,ij,nj->n', m, direction, m)[nonzero]/norms[nonzero]
    if repaired.max() > 1e-10:
        raise ArithmeticError('Finite-view separator repair failed')
    return dict(coefficients=coefficients, approximation=approximation,
        direction=direction, raw_direction=raw_direction,
        scaled_residual=float(residual), lp_objective=float(result.fun),
        dual_objective=float(dual_vector@target),
        raw_maximum_normalized_violation=float(violations.max()),
        repaired_maximum_normalized_violation=float(repaired.max()),
        direction_repair=correction, trace_signal_direction=float(np.sum(direction*s)),
        numerical_match=bool(residual < 1e-8), mixture_mass=float(coefficients.sum()),
        active_views=int(np.sum(coefficients > 1e-12)), lp_iterations=int(result.nit),
        status=int(result.status), message=result.message)


def direction_growth(direction, signal_second_moment, variance_upper=1.):
    """Feasible one-dimensional growth search, not a global optimum/upper bound.

    U=c*direction is feasible on the supplied orbit when direction is. The
    resulting T=U(I+U)^-1 has a valid moment domain for lambda_min(U)>-1/2.
    An oracle signal second moment guides this development calculation.
    """
    h = _symmetric(direction); s = _symmetric(signal_second_moment)
    v = float(variance_upper)
    if h.shape != s.shape or not np.isfinite(v) or v <= 0:
        raise ValueError('Matching matrices and positive variance upper required')
    eigenvalues, vectors = np.linalg.eigh(h)
    energy = np.diag(vectors.T@s@vectors)/v
    magnitude = np.max(np.abs(eigenvalues))
    if magnitude == 0:
        return dict(weights=h.copy(), scale=0., expected_log_lower=0., variance_upper=v)
    limit = (-.499/eigenvalues.min() if eigenvalues.min() < 0 else 1e6/magnitude)
    def objective(c):
        u = c*eigenvalues; t = u/(1+u)
        return float(energy@t+.5*np.log1p(-t*t).sum())
    scales = np.concatenate([[0.], np.geomspace(min(1e-10/magnitude, limit/100), limit, 181)])
    values = np.array([objective(c) for c in scales]); best = int(np.argmax(values))
    selected = float(scales[best]); value = float(values[best])
    if best > 0:
        fitted = minimize_scalar(lambda c: -objective(c), bounds=(scales[best-1],
            scales[min(best+1,len(scales)-1)]), method='bounded', options={'xatol': 1e-14})
        if fitted.success and -fitted.fun > value:
            selected, value = float(fitted.x), float(-fitted.fun)
    u = selected*eigenvalues; t = (vectors*(u/(1+u)))@vectors.T
    return dict(weights=t, scale=selected, expected_log_lower=max(0.,value), variance_upper=v)
