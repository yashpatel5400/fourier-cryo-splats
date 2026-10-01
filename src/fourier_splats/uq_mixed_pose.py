"""Conditional independent-centered plus common-mode pose-error bound.

Independence/centering is an assumption about errors conditional on the frozen
design, weights, density and pilot. Estimated poses do not imply it. The first
order term uses a sub-Gaussian tail, not a Gaussian approximation. All moments
are continuous cube integrals; floating-point pads are diagnostics only.
"""
import numpy as np
from .uq_cell_moments import normalized_cell_moments


def linear_moment_kernel(v):
    """Integral [1,x,y,z][1,x,y,z]' exp(2 pi i v.x) over the unit cube."""
    v = np.asarray(v, float)
    if v.shape[-1] != 3 or not np.isfinite(v).all():
        raise ValueError('Finite 3D frequencies required')
    local = normalized_cell_moments(v, 1)
    powers = np.vstack([np.zeros(3, int), np.eye(3, dtype=int)])
    out = np.empty((*v.shape[:-1], 4, 4), complex)
    for a in range(4):
        for b in range(a, 4):
            beta = powers[a]+powers[b]
            value = np.prod(np.stack([local[..., axis, beta[axis]] for axis in range(3)], axis=-1), axis=-1)
            out[..., a, b] = out[..., b, a] = value
    return out


def first_derivative_gram(k, q, coefficients, angle, shift):
    """5x5 Gram of one particle's real adjoint first derivatives, no voxel grid."""
    k, q, c = np.asarray(k, float), np.asarray(q, float), np.asarray(coefficients, complex)
    if (k.ndim != 2 or k.shape[1] != 3 or q.shape != (len(k), 2) or c.shape != (len(k),)
            or not all(np.isfinite(v).all() for v in (k, q, c))
            or not all(np.isfinite(v) and v >= 0 for v in (angle, shift))):
        raise ValueError('Finite single-particle frequencies, coefficients and scales required')
    poly = np.zeros((len(k), 5, 4), complex)
    for a in range(3):
        poly[:, a, 1:] = 2j*np.pi*angle*c[:, None]*np.cross(k, np.eye(3)[a])
    poly[:, 3:, 0] = 2j*np.pi*shift*c[:, None]*q
    minus = linear_moment_kernel(k[:, None]-k[None, :])
    plus = linear_moment_kernel(k[:, None]+k[None, :])
    gram = .5*np.real(np.einsum('qal,rbm,qrlm->ab', poly, poly.conj(), minus, optimize=True)
        +np.einsum('qal,rbm,qrlm->ab', poly, poly, plus, optimize=True))
    gram = .5*(gram+gram.T)
    # Diagnostic cancellation pad, not validated arithmetic.
    absolute = float(np.sum(abs(poly))**2)
    pad = 100*np.finfo(float).eps*max(1, len(k))*absolute
    eigen = np.linalg.eigvalsh(gram)
    if eigen[0] < -pad:
        raise ArithmeticError('First-derivative Gram is indefinite beyond diagnostic pad')
    return dict(gram=gram, spectral_norm_upper=float(np.sqrt(max(0., eigen[-1])+pad)),
        squared_norm_roundoff_diagnostic_pad=float(pad))


def mixed_pose_interval(k, q, ctf, weights, noise, angle, shift,
                        density_radius, pilot_norm, residual_norm, *,
                        alpha=.05, independent_radius=1., common_radius=0.,
                        pilot_derivative_pairings=None):
    """Return density/common/remainder bias plus a rigorous model-based tail.

    Pose vector U_i=mu_i+xi_i acts via exp(angle*[U_i[:3]]x) and
    shift*U_i[3:]. Conditional independent xi_i have E xi_i=0,
    ||xi_i||<=independent_radius_i; ||mu_i||<=common_radius_i, arbitrary
    dependence allowed for mu. xi must be independent of Gaussian measurement
    noise and of the frozen nominal design/weights. The returned bound does
    not assert any of these premises for a data-dependent estimator.
    """
    k, q, ctf = (np.asarray(v, float) for v in (k, q, ctf))
    if ctf.ndim != 2:
        raise ValueError('Particle-frequency CTF required')
    n, nq = ctf.shape
    w = np.asarray(weights, float).reshape(n, 2*nq)
    values = [noise, angle, shift, density_radius, pilot_norm, residual_norm, alpha]
    if (not np.isfinite(values).all() or noise <= 0 or not 0 < alpha < 1
            or min(values[1:6]) < 0 or k.shape != (n, nq, 3) or q.shape != (n, nq, 2)
            or not all(np.isfinite(v).all() for v in (k, q, ctf, w))):
        raise ValueError('Finite consistent design, scales and error budget required')
    independent = np.broadcast_to(independent_radius, (n,)).astype(float)
    common = np.broadcast_to(common_radius, (n,)).astype(float)
    if not all(np.isfinite(v).all() and np.all(v >= 0) for v in (independent, common)):
        raise ValueError('Finite nonnegative per-particle pose radii required')
    c = (w[:, :nq]+1j*w[:, nq:])*ctf/noise
    bounds = [first_derivative_gram(k[i], q[i], c[i], angle, shift) for i in range(n)]
    field_norms = np.array([v['spectral_norm_upper'] for v in bounds])
    R = pilot_norm+density_radius
    if pilot_derivative_pairings is None:
        scalar_bounds = R*field_norms
    else:
        pairings = np.asarray(pilot_derivative_pairings, float)
        if pairings.shape != (n, 5) or not np.isfinite(pairings).all():
            raise ValueError('Finite particle by five scaled pilot pairings required')
        scalar_bounds = np.linalg.norm(pairings, axis=1)+density_radius*field_norms
    # Uniform L2 Taylor second derivative: ||phi'^2||_2+||phi''||_2.
    kr = np.linalg.norm(k, axis=-1); qr = np.linalg.norm(q, axis=-1)
    first4 = 2*np.pi*(angle*kr*(13/720)**.25+shift*qr)
    second2 = 2*np.pi*angle**2*kr/np.sqrt(12)
    remainder_per_particle = R/2*np.sum(abs(c)*(first4**2+second2), axis=1)*(independent+common)**2
    density_bias = density_radius*residual_norm
    common_bias = float(common@scalar_bounds)
    remainder = float(remainder_per_particle.sum())
    proxy = float(np.sum((independent*scalar_bounds)**2))
    measurement = float(np.sum(w*w))
    tail = float(np.sqrt(2*np.log(2/alpha)*(measurement+proxy)))
    return dict(half_width=float(density_bias+common_bias+remainder+tail),
        density_bias=float(density_bias), common_mode_bias=common_bias,
        nonlinear_remainder_bias=remainder, stochastic_pose_variance_proxy=proxy,
        measurement_variance=measurement, subgaussian_tail=tail,
        alpha=float(alpha), independent_radius=independent.tolist(), common_radius=common.tolist(),
        scalar_derivative_bounds=scalar_bounds.tolist(), field_spectral_norm_upper=field_norms.tolist(),
        squared_norm_roundoff_diagnostic_pads=[r['squared_norm_roundoff_diagnostic_pad'] for r in bounds],
        scope='Conditional frozen-design centered independent bounded pose errors plus bounded common mode. NOT guaranteed for same-image estimated poses/weights. Ordinary floating point, no validated rounding enclosure.')
