"""Independent Gaussian noise-scale bounds with explicit calibration assumptions."""
import numpy as np
from scipy.stats import chi2


def gaussian_scale_upper(calibration, failure_probability=.01, axis=-1):
    """Bound a common SD from independent Gaussian calibration coordinates.

    Coordinates have arbitrary fixed means and common variance sigma**2. They
    must be independent of each other and of the inference noise. No fitted
    mean is subtracted: nonzero signal makes this chi-square bound conservative.
    For multiple scales, allocate failure probability across scales explicitly.
    The returned value is not valid for arbitrary colored or dependent noise.
    """
    values=np.asarray(calibration,dtype=float)
    if values.ndim==0 or not np.all(np.isfinite(values)):
        raise ValueError('Finite calibration coordinates with a sample axis required')
    if not 0<failure_probability<1:raise ValueError('Failure probability must be in (0,1)')
    df=values.shape[axis]
    if df<1:raise ValueError('At least one independent calibration coordinate required')
    return np.sqrt(np.sum(values*values,axis=axis)/chi2.ppf(failure_probability,df))
