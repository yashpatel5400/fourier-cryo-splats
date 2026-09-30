"""Fixed Gaussian row contrasts; arbitrary means and common covariance."""
import numpy as np
from scipy.linalg import helmert
from .uq_group_noise import grouped_estimator_variance_upper


def centered_grouped_variance_upper(observations, raw_weights, groups, failure_probability=.005):
    """Apply the existing noncentral bound to n-1 independent Helmert contrasts.

    The projection is deterministic. The weights must remain independent of
    the calibration/inference noises, and all prior Gaussian/common-covariance
    and exposure assumptions remain. Do not minimize this bound with another
    same-level calibration bound without accounting for multiplicity.
    """
    y = np.asarray(observations, dtype=float)
    if y.ndim != 2 or len(y) < 2 or y.shape[1] < 1 or not np.isfinite(y).all():
        raise ValueError('At least two finite calibration vectors are required')
    q = helmert(len(y), full=False)
    contrasted = q@y
    result = grouped_estimator_variance_upper(contrasted, raw_weights, groups, failure_probability)
    result.update(original_calibration_rows=len(y), independent_contrast_rows=len(y)-1,
        contrast='Deterministic orthonormal Helmert rows orthogonal to the constant vector.',
        means_assumption='Arbitrary deterministic row means; they need not be equal.')
    return result, q, contrasted
