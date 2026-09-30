"""Matrix EM with original log-kernel certificates; established mixture algebra."""
import time
import numpy as np
from scipy.special import logsumexp
from .uq_mixture_anchor import scaled_anchor_upper


def matrix_mixture_fit(log_kernels, *, weights=None, max_iterations=500,
                       tolerance=.01, check_every=25):
    """Supply finite-mixture anchors efficiently, certifying in original log units.

    Scaled matrix multiplication guides EM. Tiny entries/weights may be
    floored there, but every reported likelihood and scaled upper is evaluated
    from the original log matrix. Solver convergence is not needed for these
    real-arithmetic upper bounds. No interval-arithmetic claim is made.
    """
    begin=time.perf_counter();values=np.asarray(log_kernels,float)
    if values.ndim!=2 or min(values.shape)==0 or not np.isfinite(values).all():
        raise ValueError('Nonempty finite log-kernel matrix required')
    if max_iterations<0 or tolerance<=0 or check_every<1:raise ValueError('Invalid solver budget')
    n,cells=values.shape
    w=np.full(cells,1./cells) if weights is None else np.asarray(weights,float).copy()
    if w.shape!=(cells,) or np.any(w<0) or not np.isfinite(w).all() or w.sum()<=0:
        raise ValueError('Compatible nonnegative mixture weights required')
    w=np.maximum(w,1e-250);w/=w.sum()
    offsets=values.max(axis=1)
    kernel=np.maximum(np.exp(values-offsets[:,None]),1e-250)
    best_primal=-np.inf;best_upper=scaled_anchor_upper(values,offsets)
    best_weights=w.copy();best_upper_anchor=offsets.copy()
    history=[]
    for iteration in range(max_iterations+1):
        if iteration%check_every==0 or iteration==max_iterations:
            z=logsumexp(values+np.log(w)[None,:],axis=1)
            primal=float(z.sum());upper=scaled_anchor_upper(values,z)
            if primal>best_primal:best_primal=primal;best_weights=w.copy()
            if upper<best_upper:best_upper=upper;best_upper_anchor=z.copy()
            if best_upper<best_primal-1e-7:raise FloatingPointError('Mixture duality violation')
            gap=max(best_upper-best_primal,0.)
            history.append(dict(iteration=iteration,primal=best_primal,upper=best_upper,gap=gap))
            if gap<=tolerance or iteration==max_iterations:
                return dict(primal=best_primal,upper=best_upper,gap=gap,
                    weights=best_weights,upper_log_anchor=best_upper_anchor,
                    iterations=iteration,converged=bool(gap<=tolerance),history=history,
                    seconds=time.perf_counter()-begin)
        probability=kernel@w
        score=kernel.T@(1./probability)
        w=np.maximum(w*score,1e-250);w/=w.sum()
