"""Continuous support bounds for a total norm and negative-density norm class.

Elementary support/Jensen calculations, not an experimentally calibrated prior.
Numerical input error bounds remain the caller's responsibility.
"""
import numpy as np


def sign_class_support(positive_norm, negative_norm, radius, negative_radius):
    a, b, r, eta = map(float, [positive_norm, negative_norm, radius, negative_radius])
    if not np.isfinite([a, b, r, eta]).all() or min(a, b, r, eta) < 0 or eta > r:
        raise ValueError('Finite nonnegative norms and 0 <= negative radius <= radius required')
    total = np.hypot(a, b)
    if r*b <= eta*total:
        return float(r*total)
    return float(np.sqrt(max(0., r*r-eta*eta))*a+eta*b)


def one_sided_projection_support(norm_upper, negative_norm_lower, radius, negative_radius):
    h, b, r, eta = map(float, [norm_upper, negative_norm_lower, radius, negative_radius])
    if (not np.isfinite([h, b, r, eta]).all() or min(h, b, r, eta) < 0
            or b > h or eta > r):
        raise ValueError('Compatible norm bounds and radii required')
    return sign_class_support(np.sqrt(max(0., h*h-b*b)), b, r, eta)


def projected_sign_bounds(coefficients, norm_upper, radius, negative_radius, projection_error=0.):
    """Bounds from approximate coefficients with a supplied L2 error bound.

    The true cell-projection coefficients must be within projection_error in
    Euclidean norm. If no verified numerical error bound is supplied, the
    mathematical guarantee applies to exact coefficients/real arithmetic only.
    """
    p = np.asarray(coefficients, float); error = float(projection_error)
    if p.ndim != 1 or not np.isfinite(p).all() or not np.isfinite(error) or error < 0:
        raise ValueError('Finite coefficient vector and nonnegative projection error required')
    positive = max(0., float(np.linalg.norm(np.maximum(p, 0)))-error)
    negative = max(0., float(np.linalg.norm(np.minimum(p, 0)))-error)
    if np.hypot(positive, negative) > norm_upper*(1+1e-13)+1e-14:
        raise ValueError('Projection contradicts the supplied continuous norm bound')
    upper = one_sided_projection_support(norm_upper, negative, radius, negative_radius)
    lower_abs = one_sided_projection_support(norm_upper, positive, radius, negative_radius)
    return {'upper_support': upper, 'negative_support_magnitude': lower_abs,
            'center_offset': (upper-lower_abs)/2, 'bias_half_width': (upper+lower_abs)/2,
            'positive_norm_lower': positive, 'negative_norm_lower': negative,
            'projection_error_supplied': error}
