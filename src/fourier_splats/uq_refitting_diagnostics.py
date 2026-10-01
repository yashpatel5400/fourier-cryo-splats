"""Post hoc realized-design diagnostics, never a pose-set certificate.

True poses are required. These envelopes are exact in real arithmetic over
the declared continuous density ball at one pair of operators. Quadrature
remainders are analytic; floating-point pads remain diagnostics only.
"""
import numpy as np


def dephase_adjoint_weights(weights, q, shifts_A, field_A):
    """D* w for the real orthogonal dephasing operator D."""
    q = np.asarray(q, float)
    w = np.asarray(weights, float).reshape(len(q), 2*q.shape[1])
    shifts = np.asarray(shifts_A, float)
    if (q.ndim != 3 or q.shape[-1] != 2 or shifts.shape != (len(q), 2)
            or not np.isfinite(field_A) or field_A <= 0
            or not all(np.isfinite(v).all() for v in (q, w, shifts))):
        raise ValueError('Finite matching geometry and positive field required')
    nq = q.shape[1]
    c = (w[:, :nq]+1j*w[:, nq:])*np.exp(
        -2j*np.pi*np.einsum('nqa,na->nq', q, shifts)/field_A)
    return np.concatenate([c.real, c.imag], axis=1).ravel()


def realized_class_envelopes(true_gram, fitted_gram, weights, effective_weights,
                             centers, signs, width, density_radius, pilot_shift):
    """Worst-case total and pose-only biases at a realized operator pair.

    For T=L(rho0)+w'(D y-Ahat rho0), delta=w'(D A-Ahat)rho0,
    sup_{||h||<=B}|E T-L(rho0+h)| = |delta|+B||ell-A*D*w||.
    The expectation is conditional on the fitted operator only when inference
    noise is independent. The algebraic signal discrepancy still exists for
    same-image fitting, but does not describe that conditional expectation.
    """
    w, v = np.asarray(weights, float), np.asarray(effective_weights, float)
    B, delta = float(density_radius), float(pilot_shift)
    if (not np.isfinite([B, delta]).all() or B < 0 or w.shape != v.shape
            or not all(np.isfinite(a).all() for a in (w, v))
            or true_gram.order != fitted_gram.order
            or any(not np.array_equal(a, b) for a, b in zip(true_gram.nodes, fitted_gram.nodes))):
        raise ValueError('Matching finite weights/quadratures and nonnegative radius required')
    ft, ff = true_gram.field(v), fitted_gram.field(w)
    a, ell2 = true_gram.target(centers, signs, width)
    cross = float(v@a)
    norm2 = float(true_gram.weights@(ft*ft))
    residual2 = float(ell2-2*cross+norm2)
    pose2 = float(true_gram.weights@((ft-ff)**2))
    true_total = float(np.abs(true_gram.coefficients(v)).sum())
    fitted_total = float(np.abs(fitted_gram.coefficients(w)).sum())
    true_error = true_gram.kernel_error*true_total**2
    # A union frequency bound controls mixed terms as well as both diagonals.
    maximum = np.maximum(np.max(abs(true_gram.k), axis=(0, 1)),
                         np.max(abs(fitted_gram.k), axis=(0, 1)))
    from scipy.special import gammaln
    m = true_gram.order
    constant = 4*gammaln(m+1)-np.log(2*m+1)-3*gammaln(2*m+1)
    kernel_error = sum(np.exp(constant+2*m*np.log(4*np.pi*x)) if x else 0. for x in maximum)
    pose_error = float(kernel_error*(true_total+fitted_total)**2)
    scale = max(ell2+2*abs(cross)+norm2, (true_total+fitted_total)**2, 1.)
    pad = float(100*np.finfo(float).eps*len(w)*scale)
    if residual2 < -pad-true_error:
        raise ArithmeticError('Negative residual beyond diagnostic tolerance')
    residual = float(np.sqrt(max(0., residual2)+true_error+pad))
    pose = float(np.sqrt(pose2+pose_error+pad))
    return dict(pilot_shift=delta, true_operator_residual_norm=residual,
        operator_difference_adjoint_norm=pose,
        realized_total_bias_envelope=float(abs(delta)+B*residual),
        realized_pose_bias_envelope=float(abs(delta)+B*pose),
        residual_norm2_unpadded=residual2, pose_norm2_unpadded=pose2,
        residual_quadrature_squared_norm_bound=float(true_error),
        pose_quadrature_squared_norm_bound=pose_error,
        squared_norm_roundoff_diagnostic_pad=pad)
