"""Matched continuous generalized Gaussian-prior functional inference.

White covariance tau^2 I defines an isonormal process on L2, not an L2-valued
density. Finite observation/functional joint distributions are well defined.
These fixed-design calculations do not validate data-dependent pose estimates.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg
from scipy.stats import norm

from .uq_cell_moments import cell_fourier_moments


def cell_pose_jacobian(k, q, ctf, coefficients, box, noise):
    """Continuous-cell pilot derivatives: radians and field-fraction shifts."""
    k, q, ctf = (np.asarray(v, float) for v in (k, q, ctf))
    if (ctf.ndim != 2 or k.shape != (*ctf.shape, 3)
            or q.shape != (*ctf.shape, 2) or not np.isfinite(noise) or noise <= 0
            or not all(np.isfinite(v).all() for v in (k, q, ctf))):
        raise ValueError('Finite matching geometry and positive noise required')
    moments = cell_fourier_moments(k, coefficients, box, nthreads=1).conj()
    directions = np.stack([np.cross(k, np.eye(3)[a]) for a in range(3)], axis=-2)
    rotation = -2j*np.pi*np.einsum('nqad,nqd->nqa', directions, moments[..., 1:4])
    translation = -2j*np.pi*q*moments[..., 0, None]
    jacobian = np.concatenate([rotation, translation], axis=-1)*ctf[..., None]/noise
    return np.concatenate([jacobian.real, jacobian.imag], axis=1)


def gaussian_ball_coverage(noise_sd, bias_bound, half_width):
    """Worst fixed-bias coverage for a scalar Gaussian error and |bias|<=b."""
    s, b, h = map(float, (noise_sd, bias_bound, half_width))
    if not all(np.isfinite(v) and v >= 0 for v in (s, b, h)):
        raise ValueError('Finite nonnegative scales required')
    if s == 0:
        return float(h >= b)
    return float(norm.cdf((h-b)/s)-norm.cdf((-h-b)/s))


def continuous_gaussian(gram, centers, signs, width, tau, *, alpha=.05,
                        density_radius=None, nuisance=None, cg_rtol=1e-10,
                        maxiter=2000):
    """Solve (AA* + Sigma/tau^2)w=A ell without a density dictionary.

    Sigma=I+JJ* when a fixed pilot linearization is supplied; otherwise I.
    The returned padded variational variance is w'Sigma w+tau^2||ell-A*w||^2.
    It upper-bounds the exact posterior variance for any weights in real
    arithmetic, but is a credible interval about the posterior mean only at
    the solution. CG status/residual and the identity discrepancy are exposed.
    Quadrature has its existing real-arithmetic remainder; rounding/NUFFT do
    not have a validated enclosure. No silent fallback or outcome filtering.
    """
    tau = float(tau)
    B = tau if density_radius is None else float(density_radius)
    if (not np.isfinite(tau) or tau <= 0 or not np.isfinite(B) or B <= 0
            or not 0 < alpha < 1 or cg_rtol <= 0 or maxiter < 1):
        raise ValueError('Positive finite prior/class scales and valid solver settings required')
    a, ell2 = gram.target(centers, signs, width)
    if not np.isfinite(a).all() or not np.isfinite(ell2) or ell2 <= 0:
        raise ValueError('Finite nonzero target required')
    covariance = (lambda x: x) if nuisance is None else nuisance.covariance_apply
    ridge = tau**-2
    system = LinearOperator(gram.shape,
        matvec=lambda x: gram.matvec(x)+ridge*covariance(x), dtype=float)
    pre = gram.preconditioner(ridge) if hasattr(gram, 'preconditioner') else LinearOperator(
        gram.shape, matvec=lambda x: x/(gram.diagonal+ridge), dtype=float)
    iterations = [0]
    def count(_):
        iterations[0] += 1
    w, info = cg(system, a, M=pre, rtol=cg_rtol, atol=0., maxiter=maxiter, callback=count)
    gw = gram.matvec(w)
    error = getattr(gram, 'quadrature_error', lambda x: {'squared_field_norm': 0.,
        'gram_action_norm': 0.})(w)
    residual2 = float(ell2-2*w@a+w@gw)
    roundoff = float(50*np.finfo(float).eps*len(w)*(ell2+2*abs(w@a)+abs(w@gw)))
    if residual2 < -roundoff-error['squared_field_norm']:
        raise ArithmeticError('Negative continuous residual norm exceeds diagnostic pad')
    padded2 = max(0., residual2)+error['squared_field_norm']+roundoff
    measurement_variance = float(w@w)
    covariance_variance = float(w@covariance(w))
    pose_variance = covariance_variance-measurement_variance
    raw_variance = covariance_variance+tau**2*residual2
    padded_variance = covariance_variance+tau**2*padded2
    schur_variance = float(tau**2*(ell2-a@w))
    algebra_error = abs(raw_variance-schur_variance)/max(abs(raw_variance), 1e-300)
    solve_residual = float(np.linalg.norm(gw+ridge*covariance(w)-a)/np.linalg.norm(a))
    half = float(norm.ppf(1-alpha/2)*np.sqrt(padded_variance))
    fixed_bias = B*np.sqrt(padded2)
    return {'weights': w, 'prior_directional_sd': tau, 'density_radius': B,
        'target_norm': float(np.sqrt(ell2)), 'residual_norm': float(np.sqrt(padded2)),
        'residual_norm2_unpadded': residual2, 'measurement_variance': measurement_variance,
        'local_pose_variance': pose_variance, 'prior_residual_variance': float(tau**2*padded2),
        'variational_variance_upper': float(padded_variance),
        'schur_posterior_variance_unpadded': schur_variance,
        'variance_identity_relative_error': float(algebra_error),
        'half_width': half, 'fixed_pose_bias_bound': float(fixed_bias),
        'fixed_pose_uniform_class_coverage': gaussian_ball_coverage(np.sqrt(measurement_variance), fixed_bias, half),
        'cg_info': int(info), 'cg_iterations': iterations[0],
        'cg_relative_residual': solve_residual,
        'converged': bool(info == 0 and solve_residual <= max(2*cg_rtol, 1e-12)),
        'squared_norm_roundoff_diagnostic_pad': roundoff,
        'quadrature_errors': {k: float(v) for k, v in error.items()},
        'scope': 'Continuous fixed-design generalized Gaussian prior. Pose marginalization, if present, is pilot-linearized only. No estimated-pose or nonlinear-pose guarantee.'}
