"""Conservative Gaussian calibration allowing unknown signal means.

The common covariance and independence assumptions are essential. These helpers
do not establish those assumptions for archive particles or fitted poses.
"""
import numpy as np
from scipy.optimize import brentq


def common_covariance_trace_upper(observations, failure_probability=.005):
    """Upper bound tr(Sigma) from independent N(mu_i,Sigma) vectors.

    For t in (0,1) solving n*(t-1-log(t))/2 = log(1/beta), the bound is
    sum_i ||Y_i||^2/(n*t). Unknown means need not be equal or zero. This follows
    from the noncentral Gaussian Laplace transform and prod(1+x_j)>=1+sum x_j.
    It implies Sigma <= bound*I, usually very conservatively in many dimensions.
    """
    y=np.asarray(observations,dtype=float)
    if y.ndim!=2 or min(y.shape)<1 or not np.isfinite(y).all():
        raise ValueError('Finite nonempty matrix of independent observations required')
    if not 0<failure_probability<1:
        raise ValueError('Failure probability must lie in (0,1)')
    n=y.shape[0];target=2*np.log(1/failure_probability)/n
    # Solve for log(t) to avoid loss of significance near zero.
    log_t=brentq(lambda x:np.expm1(x)-x-target,-max(2.,2*target+2),0.)
    t=float(np.exp(log_t));energy=float(np.sum(y*y))
    return {'covariance_trace_upper':energy/(n*t),'lower_tail_fraction':t,
            'observations':n,'dimension':y.shape[1],'observed_total_energy':energy,
            'failure_probability':failure_probability,
            'assumptions':'Independent Gaussian vectors with one common covariance and arbitrary fixed means; real arithmetic.'}
