"""Fixed-estimator variance calibration allowing arbitrary within-group noise."""
import numpy as np
from .uq_noise_calibration import estimator_variance_upper


def grouped_estimator_variance_upper(observations, raw_weights, groups, failure_probability=.005):
    """Common Gaussian marginal covariance, arbitrary dependence within a group.

    Each inference group is jointly Gaussian and independent of other groups.
    Calibration rows are independent representatives from other groups. Cauchy--
    Schwarz gives Var(sum of n_g contributions) <= n_g sum of their variances.
    It is enough to calibrate the trace after multiplying weight row i by
    sqrt(n_g(i)). This can be conservative; independence within groups is not
    assumed. Weights must be independent of the calibration/inference noises.
    """
    weights = np.asarray(raw_weights, float); groups = np.asarray(groups)
    if groups.ndim != 1 or weights.ndim != 2 or len(groups) != len(weights):
        raise ValueError('One exposure label per inference weight row is required')
    labels, inverse, counts = np.unique(groups, return_inverse=True, return_counts=True)
    sizes = counts[inverse]
    result = estimator_variance_upper(observations, weights*np.sqrt(sizes[:, None]), failure_probability)
    result.update(inference_groups=len(labels), maximum_group_size=int(counts.max()),
                  inference_row_group_sizes=sizes.tolist(),
                  assumptions='Independent jointly Gaussian inference groups; arbitrary within-group dependence; independent Gaussian calibration representatives; one common marginal raw-coordinate covariance; weights fixed independently of calibration/inference noise.')
    return result


def cell_density_distance(pilot, pilot_box, reference, reference_box):
    """Exact real-arithmetic L2 distance between nested or nonnested cell bases."""
    a = np.asarray(pilot, float).reshape((pilot_box,)*3)
    b = np.asarray(reference, float).reshape((reference_box,)*3)
    ea = np.linspace(0., 1., pilot_box+1); eb = np.linspace(0., 1., reference_box+1)
    overlap = np.maximum(0., np.minimum(eb[1:, None], ea[None, 1:])-
                        np.maximum(eb[:-1, None], ea[None, :-1]))*np.sqrt(pilot_box*reference_box)
    projected = np.einsum('ia,jb,kc,abc->ijk', overlap, overlap, overlap, a, optimize=True)
    squared = np.sum(a*a)+np.sum(b*b)-2*np.sum(projected*b)
    pad = 64*np.finfo(float).eps*(np.sum(a*a)+np.sum(b*b)+2*abs(np.sum(projected*b)))
    if squared < -pad:
        raise FloatingPointError('Negative cross-cell distance')
    return float(np.sqrt(max(0., squared))), float(pad)
