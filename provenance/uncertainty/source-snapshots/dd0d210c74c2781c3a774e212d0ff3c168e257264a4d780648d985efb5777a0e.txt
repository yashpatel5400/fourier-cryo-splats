"""Integrated cube-moment bounds for the continuous nonlinear pose remainder.

These formulas refine a development certificate; they do not alter the frozen
continuous-v1 implementation. Gaussian even moments upper-bound the moments of
a linear form in independent centered uniform cube coordinates. No probability
model for the unknown density is imposed: the cube integrals are Lp norms of
known operator columns.
"""
import numpy as np


def continuous_third_derivative_bound(k, q, angle, shift):
    """L2(cube) upper bound per complex Fourier column, before CTF/whitening.

Uniformly over ||u||_2 <= 1 and t in [0,1], the phase derivative coefficient
norms are at most 2*pi*angle**j*||k|| for j=1,2,3. Detector translation adds
at most 2*pi*shift*||q|| to the first derivative. Minkowski and Holder imply
||d^3 exp(i phi)/dt^3||_2 <= ||phi'||_6^3
                                  + 3||phi'||_4 ||phi''||_4 + ||phi'''||_2.
The unit cube has volume one, variance 1/12 along each coordinate, and its even
linear-form moments are bounded by those of a Gaussian of that variance.
"""
    if angle < 0 or shift < 0:
        raise ValueError('Pose radii must be nonnegative')
    k, q = np.asarray(k), np.asarray(q)
    if k.shape[:-1] != q.shape[:-1] or k.shape[-1] != 3 or q.shape[-1] != 2:
        raise ValueError('Expected paired (...,3) volume and (...,2) detector frequencies')
    c2 = np.sqrt(1/12)
    c4 = (3/12**2)**.25
    c6 = (15/12**3)**(1/6)
    knorm = np.linalg.norm(k, axis=-1)
    qnorm = np.linalg.norm(q, axis=-1)
    translation = 2*np.pi*shift*qnorm
    first6 = 2*np.pi*angle*knorm*c6+translation
    first4 = 2*np.pi*angle*knorm*c4+translation
    second4 = 2*np.pi*angle**2*knorm*c4
    third2 = 2*np.pi*angle**3*knorm*c2
    return first6**3+3*first4*second4+third2


def integrated_cubic_remainder(k, q, ctf, weights, noise, angle, shift,
                               density_radius, pilot_norm):
    """Uniform scalar-estimator remainder over density and joint pose balls."""
    k, q, ctf = map(np.asarray, (k, q, ctf))
    if k.ndim != 3 or ctf.shape != k.shape[:2]:
        raise ValueError('Particle-frequency arrays required')
    if noise <= 0 or density_radius < 0 or pilot_norm < 0:
        raise ValueError('Invalid whitening scale or density norm bound')
    bounds = continuous_third_derivative_bound(k, q, angle, shift)
    per_particle = np.sqrt(np.sum(abs(ctf/noise)**2*bounds**2, axis=1))
    real_weights = np.asarray(weights).reshape(len(k), -1)
    if real_weights.shape[1] != 2*k.shape[1]:
        raise ValueError('Expected real and imaginary weights for each particle')
    wnorm = np.linalg.norm(real_weights, axis=1)
    operator_remainder = (density_radius+pilot_norm)*np.dot(per_particle, wnorm)/6
    count = k.shape[1]
    coefficients = (real_weights[:, :count]+1j*real_weights[:, count:])*ctf/noise
    # The weights are known: triangle inequality in Fourier frequency avoids
    # charging high-frequency operator error to weights placed at low frequency.
    # Cauchy--Schwarz shows this is no larger than the Hilbert--Schmidt variant.
    weighted_remainder = (density_radius+pilot_norm)*np.sum(abs(coefficients)*bounds)/6
    remainder = min(operator_remainder, weighted_remainder)
    return {'remainder_bias': float(remainder),
            'hilbert_schmidt_remainder_bias': float(operator_remainder),
            'frequency_weighted_remainder_bias': float(weighted_remainder),
            'per_particle_operator_third_bound': per_particle,
            'cube_linear_form_constants': {'L2': float(np.sqrt(1/12)),
                                          'L4': float((3/12**2)**.25),
                                          'L6': float((15/12**3)**(1/6))}}
