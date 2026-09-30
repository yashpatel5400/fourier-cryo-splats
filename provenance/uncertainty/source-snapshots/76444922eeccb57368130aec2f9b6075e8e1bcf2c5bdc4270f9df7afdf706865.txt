"""Fourier-cancellation Taylor remainder on an enclosing spatial ball.

Exact identities are in real arithmetic. Magnitude-based floating-point guards
are diagnostic safeguards, not validated interval bounds for special functions.
The density remains arbitrary L2 on the full unit cube.
"""
import numpy as np
from scipy.special import spherical_jn


def ball_fourier_kernel(v, radius):
    """Integral of exp(2 pi i v.x) over the centered three-dimensional ball."""
    v = np.asarray(v, dtype=float)
    if v.shape[-1] != 3 or not np.isfinite(v).all() or not np.isfinite(radius) or radius <= 0:
        raise ValueError('Finite 3D frequencies and a positive radius required')
    z = 2*np.pi*radius*np.linalg.norm(v, axis=-1)
    small = np.abs(z) < .05
    z2 = z[small]**2
    value = np.empty_like(z)
    value[small] = 1/3-z2/30+z2**2/840-z2**3/45360+z2**4/3991680
    value[~small] = spherical_jn(1, z[~small])/z[~small]
    return 4*np.pi*radius**3*value


def cube_fourier_kernel(v, radius):
    """Integral of exp(2 pi i v.x) over [-radius,radius]^3."""
    v = np.asarray(v, dtype=float)
    if v.shape[-1] != 3 or not np.isfinite(v).all() or not np.isfinite(radius) or radius <= 0:
        raise ValueError('Finite 3D frequencies and a positive radius required')
    return np.prod(2*radius*np.sinc(2*radius*v), axis=-1)


def ball_sobolev_norms(k, coefficients, radius, orders=(1, 2, 3), domain='ball'):
    """Frobenius derivative-tensor L2 norms; optional enclosing cube domain."""
    k = np.asarray(k, float); c = np.asarray(coefficients, complex)
    if k.ndim != 2 or k.shape != (len(c), 3) or not np.isfinite(k).all() or not np.isfinite(c).all():
        raise ValueError('Finite matching Fourier arrays required')
    if any(not isinstance(m, (int, np.integer)) or m < 0 for m in orders):
        raise ValueError('Nonnegative integer derivative orders required')
    if domain not in {'ball', 'cube'}:
        raise ValueError('Choose ball or cube')
    kernel = ball_fourier_kernel if domain == 'ball' else cube_fourier_kernel
    difference = kernel(k[:, None]-k[None, :], radius)
    addition = kernel(k[:, None]+k[None, :], radius)
    dot = k@k.T
    cc_conj = c[:, None]*c.conj()[None, :]
    cc = c[:, None]*c[None, :]
    norms, diagnostics = [], []
    for m in orders:
        factor = .5*(2*np.pi)**(2*m)
        terms1 = dot**m*cc_conj*difference
        terms2 = (-1)**m*dot**m*cc*addition
        squared = factor*float(np.real(np.sum(terms1)+np.sum(terms2)))
        mass = factor*float(np.sum(abs(terms1))+np.sum(abs(terms2)))
        pad = 128*np.finfo(float).eps*max(1, len(c))*mass
        if not np.isfinite(squared+pad) or squared+pad < 0:
            raise FloatingPointError('Invalid padded Sobolev norm')
        norms.append(np.sqrt(squared+pad))
        diagnostics.append({'order': m, 'squared_unpadded': squared,
                            'absolute_sum': mass, 'roundoff_pad': pad})
    return np.asarray(norms), diagnostics


def particle_ball_remainder(k, q, coefficients, angle, shift, joint_ball=True, domain='ball'):
    """Uniform L2 remainder after second-order pose expansion, one particle.

    The row-frequency pose convention is k @ exp(angle [u]_cross), with
    exp(2 pi i shift q.v) translation. Approximate detector embeddings receive
    an explicit phase-residual bound; they are never assumed exactly isometric.
    """
    k = np.asarray(k, float); q = np.asarray(q, float); c = np.asarray(coefficients, complex)
    if q.shape != (len(k), 2) or not np.isfinite(q).all():
        raise ValueError('Finite two-dimensional detector coordinates required')
    if not np.isfinite([angle, shift]).all() or angle < 0 or shift < 0:
        raise ValueError('Finite nonnegative pose radii required')
    embedding, _, rank, _ = np.linalg.lstsq(q, k, rcond=None)
    if rank != 2:
        raise ValueError('Detector coordinates must have rank two')
    singular = np.linalg.svd(embedding, compute_uv=False)
    if singular[-1] <= 1e-12*singular[0]:
        raise ValueError('Degenerate embedding')
    dual = np.linalg.pinv(embedding)
    translation_radius = shift*np.linalg.norm(dual, 2)
    cube_radius = np.sqrt(3)/2
    ball_radius = cube_radius+translation_radius
    # ||(R-I)x|| <= 2 sin(theta/2)||x|| for theta <= pi.
    # Thus the entire transformed cube lies in this expanded coordinate cube.
    expanded_cube_radius = .5+2*np.sin(min(angle, np.pi)/2)*cube_radius+translation_radius
    domain_radius = ball_radius if domain == 'ball' else expanded_cube_radius
    speed = np.hypot(angle*cube_radius, translation_radius) if joint_ball else angle*cube_radius+translation_radius
    acceleration = angle**2*cube_radius
    jerk = angle**3*cube_radius
    sobolev, diagnostics = ball_sobolev_norms(k, c, domain_radius, domain=domain)
    main = (speed**3*sobolev[2]+3*speed*acceleration*sobolev[1]+jerk*sobolev[0])/6
    # If k != q E, the dual translation produces a small extra phase.
    # Differentiate (exp(i t delta)-1) M(t) three times. |t|<=1.
    eps = 2*np.pi*shift*np.linalg.norm(q-k@dual, axis=1)
    kr = 2*np.pi*np.linalg.norm(k, axis=1)
    v, h, j = kr*speed, kr*acceleration, kr*jerk
    residual = float(np.sum(abs(c)*(eps*(v**3+3*v*h+j)+3*eps*(v*v+h)+3*eps**2*v+eps**3))/6)
    total = float(main+residual)
    if not np.isfinite(total):
        raise FloatingPointError('Nonfinite ball remainder')
    return {'field_remainder': total, 'main_remainder': float(main),
            'embedding_residual_remainder': residual,
            'maximum_phase_residual': float(np.max(eps)),
            'embedding_min_singular_value': float(singular[-1]),
            'translation_radius': float(translation_radius), 'ball_radius': float(ball_radius),
            'domain': domain, 'domain_radius': float(domain_radius),
            'path_speed_bound': float(speed), 'joint_ball': bool(joint_ball),
            'sobolev_norms': sobolev.tolist(), 'sobolev_diagnostics': diagnostics}
