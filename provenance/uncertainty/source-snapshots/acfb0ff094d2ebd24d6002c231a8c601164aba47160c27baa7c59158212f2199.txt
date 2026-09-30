"""Prototype likelihood envelopes; real-arithmetic bounds, not interval arithmetic.

No continuous rotation partition or cryo-EM noise model is validated here.
Callers must provide valid cell likelihood envelopes and an independent,
normalized predictive density before interpreting a likelihood ratio as an
e-variable. The convex optimization is standard mixture likelihood algebra.
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import logsumexp


def gaussian_ball_log_upper(residual_norm, dimension, mean_radius,
                            sigma_lower, sigma_upper):
    """Upper Gaussian density on a mean ball and a standard-deviation interval.

    A zero lower endpoint denotes a limit over positive standard deviations.
    An infinite upper endpoint is permitted. Infinite density when the ball
    contains the observation and sigma can approach zero is retained.
    """
    if not isinstance(dimension, (int, np.integer)) or dimension <= 0:
        raise ValueError('dimension must be a positive integer')
    r, b, lo, hi = np.broadcast_arrays(np.asarray(residual_norm, float),
        np.asarray(mean_radius, float), np.asarray(sigma_lower, float),
        np.asarray(sigma_upper, float))
    if (np.any(~np.isfinite(r)) or np.any(~np.isfinite(b)) or
        np.any(~np.isfinite(lo)) or np.any(np.isnan(hi)) or
        np.any(r < 0) or np.any(b < 0) or np.any(lo < 0) or
        np.any(hi <= 0) or np.any(hi < lo)):
        raise ValueError('nonnegative finite distances and valid scale interval required')
    distance = np.maximum(r-b, 0.)
    sigma = np.clip(distance / np.sqrt(dimension), lo, hi)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        value = -dimension*(np.log(sigma) + .5*np.log(2*np.pi)) - .5*(distance/sigma)**2
    return np.where(sigma == 0, np.inf, value)


@dataclass
class MixtureEnvelopeFit:
    envelope_primal: float
    log_likelihood_upper: float
    dual_gap: float
    weights: np.ndarray
    upper_anchor_weights: np.ndarray
    iterations: int
    converged: bool
    status: str
    numerical_scope: str = 'Floating-point evaluation of a real-arithmetic bound; not validated interval arithmetic.'


def mixture_envelope_upper(log_envelopes, *, tolerance=1e-6,
                           max_iterations=10000):
    """Optimize simplex weights and retain a replayable tangent upper bound.

    Each row belongs to one independent image and each column to one shared
    latent cell. The primal concerns the finite envelope relaxation; it is
    not necessarily a lower bound on the true continuous likelihood.
    """
    values = np.asarray(log_envelopes, float)
    if values.ndim != 2 or min(values.shape) == 0 or np.any(np.isnan(values)):
        raise ValueError('a nonempty image-by-cell log-envelope matrix is required')
    if not np.isfinite(tolerance) or tolerance <= 0 or max_iterations < 0:
        raise ValueError('positive tolerance and nonnegative iteration limit required')
    n, cells = values.shape
    weights = np.full(cells, 1./cells)
    if np.any(np.isposinf(values)):
        return MixtureEnvelopeFit(np.inf, np.inf, np.inf, weights, weights.copy(),
                                  0, False, 'infinite_cell_envelope')
    if np.any(np.all(np.isneginf(values), axis=1)):
        raise ValueError('every image must have at least one positive envelope')
    offsets = values.max(axis=1)
    scaled = values-offsets[:, None]
    best_primal, best_upper = -np.inf, np.inf
    best_weights, upper_weights = weights.copy(), weights.copy()
    for iteration in range(max_iterations+1):
        log_weights = np.log(weights)
        log_probability = logsumexp(scaled+log_weights[None, :], axis=1)
        # Evaluate the supporting upper bound without dropping small kernels.
        log_scores = logsumexp(scaled-log_probability[:, None], axis=0)
        primal = float(np.sum(log_probability)+np.sum(offsets))
        with np.errstate(over='ignore'):
            gap = float(np.exp(np.max(log_scores))-n)
        if gap < -1e-8*max(n, 1):
            raise FloatingPointError('Numerical mixture duality violation')
        upper = primal + max(gap, 0.)
        if primal > best_primal:
            best_primal, best_weights = primal, weights.copy()
        if upper < best_upper:
            best_upper, upper_weights = upper, weights.copy()
        remaining_gap = max(best_upper-best_primal, 0.)
        if remaining_gap <= tolerance or iteration == max_iterations:
            return MixtureEnvelopeFit(best_primal, best_upper, remaining_gap,
                best_weights, upper_weights, iteration, remaining_gap <= tolerance,
                'converged' if remaining_gap <= tolerance else 'iteration_limit')
        # EM for the finite nonnegative likelihood matrix. Tiny positive weights
        # avoid silently removing columns through floating-point underflow.
        log_weights += log_scores
        log_weights -= logsumexp(log_weights)
        weights = np.maximum(np.exp(log_weights), np.finfo(float).tiny)
        weights /= weights.sum()


def log_validation_evalue(log_predictive_density, log_likelihood_upper):
    """Arithmetic only: normalization/independence are obligations of the caller."""
    q, upper = float(log_predictive_density), float(log_likelihood_upper)
    if not np.isfinite(q) or np.isnan(upper) or np.isneginf(upper):
        raise ValueError('finite predictive log density and valid upper bound required')
    return q-upper
