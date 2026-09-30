"""Choose a feasible convex mixture of spectral dual support directions."""
import numpy as np
from scipy.optimize import minimize


def closest_support_mixture(directions, target, initial=None):
    """Minimize distance to a convex hull; every returned mixture is feasible.

    This small quadratic program tightens a dual residual, not the primal
    interval. Solver success is not required for validity. A projected feasible
    candidate is retained only when no worse than the supplied initial mixture.
    """
    directions, target = np.asarray(directions, float), np.asarray(target, float)
    if directions.ndim != 2 or target.shape != (directions.shape[1],) or len(directions) == 0:
        raise ValueError('Expected nonempty direction rows and a compatible target')
    if not np.isfinite(directions).all() or not np.isfinite(target).all():
        raise ValueError('Finite directions and target required')
    count = len(directions)
    initial = np.full(count, 1/count) if initial is None else np.asarray(initial, float).copy()
    if initial.shape != (count,) or not np.isfinite(initial).all() or (initial < 0).any() or initial.sum() <= 0:
        raise ValueError('Nonnegative nonzero initial mixture required')
    initial /= initial.sum()
    gram = directions@directions.T; linear = directions@target
    scale = max(float(np.max(abs(gram))), float(np.max(abs(linear))), np.finfo(float).tiny)
    if count == 1:
        candidate = initial; success = True; message = 'One available support direction'
    else:
        result = minimize(lambda x: (.5*x@gram@x-linear@x)/scale, initial,
            jac=lambda x: (gram@x-linear)/scale, method='SLSQP', bounds=[(0., 1.)]*count,
            constraints=[{'type': 'eq', 'fun': lambda x: x.sum()-1., 'jac': lambda x: np.ones_like(x)}],
            options={'ftol': 1e-13, 'maxiter': 500})
        candidate = np.maximum(result.x, 0.)
        candidate = candidate/candidate.sum() if np.isfinite(candidate).all() and candidate.sum() > 0 else initial
        success, message = bool(result.success), str(result.message)
    old = float(np.linalg.norm(target-initial@directions))
    new = float(np.linalg.norm(target-candidate@directions))
    if new > old: candidate, new = initial, old
    return {'mixture': candidate, 'initial_distance': old, 'distance': new,
            'optimizer_success': success, 'optimizer_message': message}
