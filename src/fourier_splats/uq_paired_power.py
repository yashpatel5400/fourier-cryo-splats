"""Bilinear Gaussian ingredients for a proposed paired-exposure map test.

These routines do not certify a map's continuous orientation orbit, the noise
covariance bound, or exposure independence. They are not a complete uncertainty
method. A valid e-value needs the external null-exponent bound in the proposal.
"""
import numpy as np


def _parameters(weights, variance_upper):
    t, v = np.broadcast_arrays(np.asarray(weights, float), np.asarray(variance_upper, float))
    if not np.isfinite(t).all() or not np.isfinite(v).all() or np.any(v < 0):
        raise ValueError('Finite weights and nonnegative variance bounds required')
    if np.any(np.abs(t)*v >= 1):
        raise ValueError('Bilinear Gaussian moment requires abs(weight)*variance < 1')
    return t, v


def complex_log_normalizer(weights, variance_upper):
    """Additive log factor cancelling the bounding zero-mean Gaussian moment.

    Each frequency has two real coordinates, each of variance at most v in
    Loewner order. This is not the convention E[abs(complex noise)**2] = v.
    """
    t, v = _parameters(weights, variance_upper)
    return np.log1p(-(t*v)**2)


def common_mean_coefficient(weights, variance_upper):
    """Coefficient of squared complex mean at the bounding covariance."""
    t, v = _parameters(weights, variance_upper)
    return t/(1-t*v)


def mismatch_coefficient(weights, variance_upper, gain_low=1., gain_high=1., phase_radius=0.):
    """Upper coefficient with mean Z = gain * exp(i phase) * mean X.

    Gains are nonnegative in [low, high]. The phase belongs to [-radius,
    radius]. All mismatches may differ by frequency; their separate maxima
    form a sufficient upper bound even if physical parameters couple them.
    """
    t, v = _parameters(weights, variance_upper)
    t,v,lo,hi,rad = np.broadcast_arrays(t,v,np.asarray(gain_low,float),
        np.asarray(gain_high,float),np.asarray(phase_radius,float))
    if not all(np.isfinite(x).all() for x in [lo,hi,rad]) or np.any(lo<0) or np.any(hi<lo) or np.any(rad<0):
        raise ValueError('Finite ordered nonnegative gain bounds and phase radii required')
    cosine = np.where(t>=0,1.,np.cos(np.minimum(rad,np.pi)))
    denominator=1-(t*v)**2
    def coefficient(gain):
        return (t*gain*cosine+.5*t*t*v*(1+gain*gain))/denominator
    # The numerator is convex in gain, so its maximum occurs at an endpoint.
    return np.maximum(coefficient(lo),coefficient(hi))


def paired_power_log_factor(first, second, weights, variance_upper, null_exponent_upper):
    """Evaluate a factor given an externally certified noncentral bound.

    The supplied scalar bound must dominate the noncentral exponent under
    every null map/pose/amplitude. Zero is valid for a certified power-cone
    constraint, not for an arbitrary candidate. Parameters must be fixed
    independently of these inference images (or predictably across particles).
    """
    x,z,t,v=np.broadcast_arrays(np.asarray(first,complex),np.asarray(second,complex),
        np.asarray(weights,float),np.asarray(variance_upper,float))
    if x.ndim!=1 or not np.isfinite(x).all() or not np.isfinite(z).all():
        raise ValueError('One finite complex frequency vector per exposure required')
    if not np.isfinite(null_exponent_upper) or null_exponent_upper<0:
        raise ValueError('An externally certified finite nonnegative exponent bound is required')
    return float(np.sum(t*np.real(x*np.conj(z))+complex_log_normalizer(t,v))-null_exponent_upper)
