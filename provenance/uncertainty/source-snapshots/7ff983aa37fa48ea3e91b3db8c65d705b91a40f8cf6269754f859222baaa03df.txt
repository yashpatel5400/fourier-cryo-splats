"""Shared Gaussian scale via a convex-in-precision mixture term.

Real-arithmetic bounds evaluated in ordinary floating point. This prototype
 does not validate continuous poses or the experimental Gaussian noise model.
"""
import numpy as np
from .uq_mixture_validation import mixture_envelope_upper


def profile_common_scale(residual_norms, dimension, *, tolerance=1.,
                         max_nodes=65, mixture_iterations=2000,
                         mixture_tolerance=1e-4):
    r = np.asarray(residual_norms, float)
    if r.ndim != 2 or min(r.shape) == 0 or not np.isfinite(r).all() or np.any(r < 0):
        raise ValueError('Finite nonnegative image-by-cell residual distances required')
    if not isinstance(dimension, (int, np.integer)) or dimension < 1:
        raise ValueError('Positive integer dimension required')
    if tolerance <= 0 or not np.isfinite(tolerance) or max_nodes < 1:
        raise ValueError('Positive gap target and node budget required')
    lower_scale, upper_scale = float(r.min()/np.sqrt(dimension)), float(r.max()/np.sqrt(dimension))
    if lower_scale == 0:
        return dict(log_likelihood_upper=np.inf, feasible_log_likelihood=-np.inf,
            gap=np.inf, converged=False, status='zero_distance_unresolved',
            bracket=[lower_scale, upper_scale], nodes=[], leaves=[], points=[])
    coefficient = .5*len(r)*dimension
    constant = -coefficient*np.log(2*np.pi)
    squared = r*r
    points, cache, nodes = [], {}, []
    best_lower = -np.inf
    def point(precision):
        nonlocal best_lower
        precision = float(precision)
        if precision in cache:
            return cache[precision]
        fit = mixture_envelope_upper(-.5*precision*squared,
            tolerance=mixture_tolerance, max_iterations=mixture_iterations)
        normalization = coefficient*np.log(precision)+constant
        identifier = len(points)
        points.append(dict(id=identifier, precision=precision, scale=precision**-.5,
            mixture_term_lower=fit.envelope_primal,
            mixture_term_upper=fit.log_likelihood_upper,
            feasible_log_likelihood=float(fit.envelope_primal+normalization),
            optimizer_gap=fit.dual_gap, optimizer_iterations=fit.iterations,
            optimizer_converged=fit.converged, weights=fit.weights.tolist(),
            upper_anchor_weights=fit.upper_anchor_weights.tolist()))
        best_lower = max(best_lower, points[-1]['feasible_log_likelihood'])
        cache[precision] = identifier
        return identifier
    initial_scale = float(np.clip(1., lower_scale, upper_scale))
    initial_id = point(initial_scale**-2)
    def evaluate(lo, hi, parent):
        # The supremum over weights of the convex log-sum-exp term is convex
        # in precision. Endpoint upper bounds give a global chord upper bound.
        a, b = float(hi**-2), float(lo**-2)
        left, right = point(a), point(b)
        fa, fb = points[left]['mixture_term_upper'], points[right]['mixture_term_upper']
        if a == b:
            maximizing = a
            upper = fa+coefficient*np.log(a)+constant
        else:
            slope = (fb-fa)/(b-a)
            maximizing = float(np.clip(-coefficient/slope, a, b)) if slope < 0 else b
            upper = fa+slope*(maximizing-a)+coefficient*np.log(maximizing)+constant
        candidate = point(maximizing)
        identifier = len(nodes)
        nodes.append(dict(id=identifier, parent=parent, scale_interval=[lo,hi],
            precision_interval=[a,b], endpoint_point_ids=[left,right],
            bound_maximizer_point_id=candidate, log_likelihood_upper=float(upper)))
        return identifier
    leaves = [evaluate(lower_scale, upper_scale, None)]
    while True:
        largest = max(leaves, key=lambda i: nodes[i]['log_likelihood_upper'])
        upper = nodes[largest]['log_likelihood_upper']
        gap = upper-best_lower
        if gap < -1e-7:
            raise FloatingPointError('Scale upper bound below a feasible likelihood')
        finished = bool(gap <= tolerance)
        exhausted = len(nodes)+2 > max_nodes or lower_scale == upper_scale
        if finished or exhausted:
            return dict(log_likelihood_upper=upper, feasible_log_likelihood=best_lower,
                gap=max(gap,0.), converged=finished,
                status='converged' if finished else 'node_limit',
                bracket=[lower_scale,upper_scale], nodes=nodes, leaves=leaves,
                points=points, initial_point_id=initial_id,
                numerical_scope='Ordinary floating point; not validated interval arithmetic.')
        parent = nodes[largest]; lo, hi = parent['scale_interval']
        mid = float(np.exp(.5*(np.log(lo)+np.log(hi))))
        if not lo < mid < hi:
            return dict(log_likelihood_upper=upper, feasible_log_likelihood=best_lower,
                gap=max(gap,0.), converged=False, status='floating_scale_resolution_limit',
                bracket=[lower_scale,upper_scale], nodes=nodes, leaves=leaves,
                points=points, initial_point_id=initial_id)
        leaves.remove(largest)
        leaves.extend([evaluate(lo,mid,largest), evaluate(mid,hi,largest)])
