"""Revision interval numerics; frozen experiments retain their original module."""
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm


def bias_aware_half_width_stable(sd, bias, alpha=.05):
    """Folded-normal critical value, solving both tails without a bias shortcut.

    Solve sf(t)+sf(t+2*bias/sd)=alpha for t=(q-bias)/sd in log-tail space.
    An outward numerical pad covers root-solver tolerance to working precision.
    This remains floating-point numerical evaluation, not interval arithmetic.
    """
    sd, bias, alpha = float(sd), float(bias), float(alpha)
    if not np.isfinite([sd, bias, alpha]).all() or sd < 0 or bias < 0 or not 0 < alpha < 1:
        raise ValueError('Finite nonnegative scales and alpha in (0,1) required')
    if sd == 0:
        return bias
    with np.errstate(over='ignore'):
        ratio = np.float64(bias)/sd
        twice = 2*ratio
    lower = max(-ratio, norm.isf(alpha)); upper = norm.isf(alpha/2)
    log_alpha = np.log(alpha)
    def tail(t):
        return np.logaddexp(norm.logsf(t), norm.logsf(t+twice))-log_alpha
    if tail(lower) <= 0:
        root = lower  # Remote tail can be below floating-point resolution.
    elif tail(upper) >= 0:
        root = upper
    else:
        root = brentq(tail, lower, upper, xtol=np.nextafter(0., 1.), rtol=4*np.finfo(float).eps)
    q = bias+sd*(root+8*np.finfo(float).eps*max(1., abs(root)))
    return float(np.nextafter(max(0., q), np.inf))
