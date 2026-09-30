"""Selected-block continuous trace diagnostics for polynomial Fourier fields.

These analytic integration identities check scale efficiency. Ordinary floating
point here is not a validated upper bound on all operator roundoff errors.
"""
import numpy as np
from numpy.polynomial.legendre import poly2leg
from scipy.special import spherical_jn


def interval_polynomial_moments(frequency, degree=6, radius=.5):
    """Integral x^j exp(2 pi i frequency x) on [-radius,radius]."""
    frequency = np.asarray(frequency, float)
    if degree < 0 or not isinstance(degree, int) or radius <= 0 or not np.isfinite(frequency).all():
        raise ValueError('Finite frequencies, nonnegative degree and positive radius required')
    t = 2*np.pi*radius*frequency
    bessel = np.stack([spherical_jn(j, abs(t))*np.where(t < 0, (-1)**j, 1)*(1j**j) for j in range(degree+1)], axis=-1)
    out = np.empty(frequency.shape+(degree+1,), complex)
    for j in range(degree+1):
        coefficients = poly2leg(np.r_[np.zeros(j), 1.])
        out[..., j] = 2*radius*radius**j*np.sum(bessel[..., :j+1]*coefficients, axis=-1)
    return out


def polynomial_fourier_block_trace(k, coefficients, monomials):
    """Sum_a ||Re sum_qp C[q,a,p] x^beta_p exp(2pi i k_q.x)||^2.

    Uses pair frequencies and exact polynomial cube moments, rather than the
    low-order grid that chose the scales. Returns a heuristic arithmetic guard.
    """
    k = np.asarray(k, float); c = np.asarray(coefficients, complex)
    if k.ndim != 2 or k.shape[1] != 3 or c.ndim != 3 or len(c) != len(k) or c.shape[2] != len(monomials):
        raise ValueError('Matching frequencies and polynomial field coefficients required')
    if not np.isfinite(k).all() or not np.isfinite(c).all(): raise ValueError('Finite inputs required')
    max_degree = 2*max(sum(beta) for beta in monomials)
    minus = interval_polynomial_moments(k[:, None]-k[None, :], degree=max_degree)
    plus = interval_polynomial_moments(k[:, None]+k[None, :], degree=max_degree)
    active = [p for p in range(c.shape[2]) if np.any(c[:, :, p] != 0)]
    total = 0.; absolute_sum = 0.
    for pos, p in enumerate(active):
        for s in active[pos:]:
            beta = tuple(a+b for a, b in zip(monomials[p], monomials[s]))
            kernel_minus = np.prod(np.stack([minus[..., axis, power] for axis, power in enumerate(beta)]), axis=0)
            kernel_plus = np.prod(np.stack([plus[..., axis, power] for axis, power in enumerate(beta)]), axis=0)
            left = c[:, :, p]; right = c[:, :, s]
            tminus = (left@right.conj().T)*kernel_minus
            tplus = (left@right.T)*kernel_plus
            factor = .5 if p == s else 1.
            total += factor*float(np.real(np.sum(tminus)+np.sum(tplus)))
            absolute_sum += factor*float(np.sum(abs(tminus))+np.sum(abs(tplus)))
    guard = 128*np.finfo(float).eps*max(1, len(k))*max(1, len(active))*absolute_sum
    if not np.isfinite(total+guard) or total+guard < 0:
        raise FloatingPointError('Invalid continuous block trace')
    return {'trace_unpadded': total, 'heuristic_roundoff_guard': guard,
            'norm_unpadded': float(np.sqrt(max(0, total))),
            'absolute_term_sum': absolute_sum, 'active_spatial_monomials': len(active)}
