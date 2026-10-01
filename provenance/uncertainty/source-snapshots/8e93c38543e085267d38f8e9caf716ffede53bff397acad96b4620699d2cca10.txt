"""Solvable alignment-bias control, not a cryo-EM reconstruction estimator."""
import numpy as np
from scipy.optimize import brentq


def cross_power_amplitude_interval(mean_cross, particles, sigma, alpha=.05):
    """Invert a product-normal sub-gamma bound; endpoints concern amplitude."""
    if (not isinstance(particles, (int, np.integer)) or particles < 1
            or not np.isfinite([mean_cross, sigma, alpha]).all()
            or sigma <= 0 or not 0 < alpha < 1):
        raise ValueError('Finite mean, positive noise and integer sample size required')
    variance = sigma**2
    b = variance*np.log(2/alpha)/particles
    def margin(m):
        return abs(mean_cross-m)-2*np.sqrt(b*(variance+m))-b
    # The margin is convex; on its increasing-side branch its derivative
    # vanishes at m=b-variance if that point is to the right of the kink.
    mode = max(0., mean_cross, b-variance)
    if margin(mode) > 0:
        return dict(empty=True, lower=None, upper=None)
    lower = 0. if margin(0.) <= 0 else brentq(margin, 0., mode, xtol=1e-12)
    hi = max(1., mode+1.)
    while margin(hi) <= 0:
        hi *= 2
    upper = brentq(margin, mode, hi, xtol=1e-12)
    return dict(empty=False, lower=float(np.sqrt(lower)), upper=float(np.sqrt(upper)))
