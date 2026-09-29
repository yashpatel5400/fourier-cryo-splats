"""Gaussian-start power upper bounds with an explicit failure probability.

Classical randomized norm estimation, not a new general statistical theorem.
The matrix must be fixed before the fresh Gaussian probes are drawn. Numerical
roundoff/operator approximation is not certified by this probability bound.
"""
import numpy as np
from scipy.special import erfinv


def gaussian_power_upper(matvec, dimension, failure_probability=1e-6,
                         probes=4, iterations=40, seed=0, trace_upper=np.inf,
                         callback=None):
    """Upper-bound lambda_max of a fixed positive semidefinite operator.

    For independent g_j~N(0,I), c=sqrt(2)*erfinv(delta**(1/m)). With probability
    at least 1-delta, max_j |<v_max,g_j>| >= c. Consequently, simultaneously
    for every integer p>=1, lambda_max <= (max_j ||C**p g_j|| / c)**(1/p).
    A Rayleigh value is only a lower bound. Normalization prevents overflow.
    Any deterministic trace upper bound may further cap the reported upper.
    """
    if dimension < 1 or probes < 1 or iterations < 1 or not 0 < failure_probability < 1:
        raise ValueError('Positive dimensions/iterations and delta in (0,1) required')
    rng = np.random.default_rng(seed)
    vectors = rng.normal(size=(probes, dimension))
    initial_norm = np.linalg.norm(vectors, axis=1)
    vectors /= initial_norm[:, None]
    log_norm = np.log(initial_norm)
    threshold = float(np.sqrt(2)*erfinv(failure_probability**(1/probes)))
    history = []; best_upper = float(trace_upper); lower = 0.
    for step in range(1, iterations+1):
        for j in range(probes):
            image = np.asarray(matvec(vectors[j]), dtype=float)
            length = float(np.linalg.norm(image))
            lower = max(lower, float(vectors[j]@image))
            if not np.isfinite(length) or length <= 0:
                raise FloatingPointError('Degenerate/nonfinite power iterate; use deterministic audit instead')
            log_norm[j] += np.log(length)
            vectors[j] = image/length
        upper = float(np.exp((np.max(log_norm)-np.log(threshold))/step))
        best_upper = min(best_upper, upper)
        record = {'iteration': step, 'power_upper': upper, 'best_upper': best_upper,
                  'rayleigh_lower': lower, 'upper_over_lower': best_upper/lower if lower>0 else None}
        history.append(record)
        if callback is not None:
            callback(record)
    return {'eigenvalue_upper': best_upper, 'rayleigh_lower': lower,
            'failure_probability': failure_probability, 'probes': probes,
            'iterations': iterations, 'seed': seed, 'projection_threshold': threshold,
            'history': history,
            'scope': 'Probability over fresh design-only Gaussian probes for a fixed PSD matrix, in real arithmetic; not validated floating point.'}
