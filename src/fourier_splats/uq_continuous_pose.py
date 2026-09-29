"""Continuous-density nonlinear-pose sensitivity bounds on the unit cube.

The density remains arbitrary L2, not a quadrature-node vector. Polynomial
Fourier adjoint derivatives are integrated with an explicit Gauss remainder.
Spectral relaxations retain the shared density field; the initial bound uses
||rho0||+B for its pose contribution and can be conservative. All numerical
linear algebra uses ordinary floating point, not validated interval arithmetic.
"""
import sys
import numpy as np
import finufft
from numpy.polynomial.legendre import leggauss
from scipy.special import gammaln, logsumexp
from scipy.linalg import eigvalsh
from scipy.spatial.transform import Rotation
from .uq_continuous import cell_forward

PAIRS = [(a, b) for a in range(5) for b in range(a, 5)]
PAIR_SCALE = np.array([1. if a == b else np.sqrt(2.) for a, b in PAIRS])


def cube_quadrature(order):
    nodes, weights = leggauss(order)
    nodes = nodes / 2
    weights = weights / 2
    z, y, x = np.meshgrid(nodes, nodes, nodes, indexing='ij')
    xyz = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    quadrature = (weights[:, None, None] * weights[None, :, None] * weights[None, None, :]).ravel()
    return xyz, quadrature


def polynomial_kernel_error(k, order, degree=4):
    """Uniform tensor Gauss error for x**beta exp(2 pi i v.x).

    Coordinate frequencies satisfy |v_a| <= 2 max|k_a|; each beta_a <= degree.
    Positive quadrature weights and the Peano kernel justify the complex-valued
    derivative bound without adding a spurious independent real/imag factor.
    """
    if order < 1 or degree < 0:
        raise ValueError('Positive quadrature order and nonnegative degree required')
    n = 2 * order
    log_constant = 4*gammaln(order+1)-np.log(2*order+1)-3*gammaln(2*order+1)
    maximum = np.max(np.abs(k).reshape(-1, 3), axis=0)
    errors = []
    for maximum_axis in maximum:
        frequency = 4*np.pi*maximum_axis
        beta_bounds = []
        for beta in range(degree+1):
            terms = []
            for derivative in range(min(beta, n)+1):
                power = n-derivative
                if frequency == 0 and power:
                    continue
                term = (gammaln(n+1)-gammaln(derivative+1)-gammaln(n-derivative+1)
                        +gammaln(beta+1)-gammaln(beta-derivative+1)
                        +(beta-derivative)*np.log(.5))
                if power:
                    term += power*np.log(frequency)
                terms.append(term)
            beta_bounds.append(np.exp(log_constant+logsumexp(terms)) if terms else 0.)
        errors.append(max(beta_bounds))
    return float(sum(errors))


def pose_derivative_fields(k, q, coefficients, xyz, angle, shift, backend='nufft'):
    """First and symmetric second derivatives of one particle's adjoint field.

    k and q use cycles per physical field; shift uses fractions of that field.
    The rotation convention is row k @ exp(angle [u]_cross). Hessian columns
    use orthonormal symmetric vectorization: off-diagonal entries times sqrt(2).
    """
    k, q, xyz = map(np.asarray, (k, q, xyz))
    frequencies = np.concatenate([k, q], axis=1)
    moments = np.concatenate([frequencies.T, np.stack([frequencies[:, a]*frequencies[:, b] for a, b in PAIRS])])
    coefficients = np.ascontiguousarray(moments*np.asarray(coefficients)[None], dtype=complex)
    if backend == 'direct':
        transform = np.empty((len(xyz), 20), dtype=complex)
        for start in range(0, len(xyz), 4096):
            phase = np.exp(2j*np.pi*(k@xyz[start:start+4096].T))
            transform[start:start+len(phase.T)] = (coefficients@phase).T
    elif backend == 'nufft':
        source = [np.ascontiguousarray(2*np.pi*k[:, a]) for a in range(3)]
        dest = [np.ascontiguousarray(xyz[:, a]) for a in range(3)]
        transform = finufft.nufft3d3(*source, coefficients, *dest, isign=1, eps=1e-12,
                                    nthreads=1 if sys.platform == 'darwin' else 4).T
    else:
        raise ValueError('Choose direct or nufft')
    first = transform[:, :5]
    second = np.empty((len(xyz), 5, 5), complex)
    for index, (a, b) in enumerate(PAIRS):
        second[:, a, b] = second[:, b, a] = transform[:, 5+index]
    linear = np.zeros((len(xyz), 5, 5))
    x, y, z = xyz.T
    scale = 2*np.pi
    linear[:, 0, 1] = -scale*angle*z; linear[:, 0, 2] = scale*angle*y
    linear[:, 1, 0] = scale*angle*z; linear[:, 1, 2] = -scale*angle*x
    linear[:, 2, 0] = -scale*angle*y; linear[:, 2, 1] = scale*angle*x
    linear[:, 3, 3] = scale*shift; linear[:, 4, 4] = scale*shift
    d1 = (1j*np.einsum('par,pr->pa', linear, first)).real
    d2 = -np.einsum('par,prs,pbs->pab', linear, second, linear).real
    dot = np.sum(xyz*first[:, :3], axis=1)
    for a in range(3):
        for b in range(3):
            curvature = scale*angle**2*(.5*(xyz[:, b]*first[:, a]+xyz[:, a]*first[:, b])-(a == b)*dot)
            d2[:, a, b] += (1j*curvature).real
    return d1, np.stack([d2[:, a, b]*scale for (a, b), scale in zip(PAIRS, PAIR_SCALE)], axis=1)


def derivative_coefficient_bounds(k, q, coefficients, angle, shift):
    """Absolute coefficient sums for polynomial Fourier derivative columns."""
    k, q = np.asarray(k), np.asarray(q)
    first = np.zeros((len(k), 5))
    for a in range(3):
        first[:, a] = 2*np.pi*angle*(np.sum(abs(k), axis=1)-abs(k[:, a]))
    first[:, 3:] = 2*np.pi*shift*abs(q)
    second = []
    for (a, b), pair_scale in zip(PAIRS, PAIR_SCALE):
        curvature = np.zeros_like(k)
        if a < 3 and b < 3:
            curvature[:, b] += .5*k[:, a]
            curvature[:, a] += .5*k[:, b]
            if a == b:
                curvature -= k
            curvature *= 2*np.pi*angle**2
        second.append(pair_scale*(first[:, a]*first[:, b]+np.sum(abs(curvature), axis=1)))
    absolute = abs(coefficients)
    return absolute@first, absolute@np.stack(second, axis=1)


def grouped_gram_bounds(gram, coefficient_sums, block_sizes, kernel_error):
    """Triangle and shared-field spectral bounds, with integration error pads."""
    gram = np.asarray(gram); coefficient_sums = np.asarray(coefficient_sums)
    if sum(block_sizes) != len(gram):
        raise ValueError('Block sizes do not partition Gram')
    offset = np.cumsum([0]+list(block_sizes))
    squared = np.array([max(0, float(np.trace(gram[a:b, a:b]))) for a, b in zip(offset[:-1], offset[1:])])
    masses = np.array([coefficient_sums[a:b]@coefficient_sums[a:b] for a, b in zip(offset[:-1], offset[1:])])
    norms = np.sqrt(squared+kernel_error*masses)
    if not np.any(norms > 0):
        return {'triangle': 0., 'spectral': 0., 'minimum': 0., 'quadrature_eigenvalue_pad': 0.}
    # Any positive block scales are valid; an exactly zero field contributes none.
    active_blocks = norms > 0
    active_columns = np.repeat(active_blocks, block_sizes)
    scales = np.repeat(np.sqrt(norms[active_blocks]), np.array(block_sizes)[active_blocks])
    scaled = gram[np.ix_(active_columns, active_columns)]/scales[:, None]/scales[None, :]
    sums = coefficient_sums[active_columns]/scales
    integration_pad = kernel_error*float(sums@sums)
    roundoff = 20*np.finfo(float).eps*max(1, len(scaled))*np.linalg.norm(scaled, ord=np.inf)
    eigenvalue = float(eigvalsh(scaled, subset_by_index=[len(scaled)-1, len(scaled)-1], check_finite=False)[0])
    spectral = float(np.sqrt(np.sum(norms)*max(0, eigenvalue+integration_pad+roundoff)))
    triangle = float(sum(norms))
    return {'triangle': triangle, 'spectral': spectral, 'minimum': min(triangle, spectral),
            'quadrature_eigenvalue_pad': integration_pad, 'roundoff_eigenvalue_pad': float(roundoff)}


def continuous_pose_audit(k, q, ctf, weights, noise_std, angle, shift, density_radius,
                          pilot_norm, residual_norm, order=32, backend='nufft', callback=None):
    """Bound the total scalar bias for fixed, inference-noise-independent weights."""
    if not 0 <= angle <= np.pi or min(shift, density_radius, pilot_norm, residual_norm) < 0:
        raise ValueError('Invalid pose or density radii')
    k, q, ctf = map(np.asarray, (k, q, ctf)); n, nq = ctf.shape
    w = np.asarray(weights).reshape(n, 2*nq)
    coefficients = (w[:, :nq]+1j*w[:, nq:])*ctf/noise_std
    xyz, quadrature = cube_quadrature(order); sqrt_quad = np.sqrt(quadrature)
    # Columns are first-order blocks, followed by half-scaled Hessian blocks.
    fields = np.empty((len(xyz), 20*n)); sums = np.empty(20*n)
    for i in range(n):
        first, second = pose_derivative_fields(k[i], q[i], coefficients[i], xyz, angle, shift, backend)
        m1, m2 = derivative_coefficient_bounds(k[i], q[i], coefficients[i], angle, shift)
        fields[:, 5*i:5*(i+1)] = sqrt_quad[:, None]*first
        fields[:, 5*n+15*i:5*n+15*(i+1)] = .5*sqrt_quad[:, None]*second
        sums[5*i:5*(i+1)] = m1; sums[5*n+15*i:5*n+15*(i+1)] = .5*m2
        if callback is not None and ((i+1) % 16 == 0 or i+1 == n):
            callback({'particle_derivatives_completed': i+1, 'total': n})
    gram = fields.T@fields
    storage = fields.nbytes+gram.nbytes+xyz.nbytes+quadrature.nbytes+sums.nbytes
    del fields
    error = polynomial_kernel_error(k, order)
    first = grouped_gram_bounds(gram[:5*n, :5*n], sums[:5*n], [5]*n, error)
    second = grouped_gram_bounds(gram[5*n:, 5*n:], sums[5*n:], [15]*n, error)
    joint = grouped_gram_bounds(gram, sums, [5]*n+[15]*n, error)
    field_bound = min(first['minimum']+second['minimum'], joint['minimum'])
    kr = np.linalg.norm(k, axis=-1); radius = np.sqrt(3)/2
    speed = 2*np.pi*(angle*kr*radius+shift*np.linalg.norm(q, axis=-1))
    acceleration = 2*np.pi*angle**2*kr*radius
    jerk = 2*np.pi*angle**3*kr*radius
    third_norm = np.linalg.norm(abs(ctf/noise_std)*(speed**3+3*speed*acceleration+jerk), axis=1)
    remainder = float((pilot_norm+density_radius)/6*np.sum(third_norm*np.linalg.norm(w, axis=1)))
    density_bias = float(density_radius*residual_norm)
    pose_bias = float((pilot_norm+density_radius)*field_bound)
    return {'density_bias': density_bias, 'pose_polynomial_bias': pose_bias, 'remainder_bias': remainder,
            'total_bias': density_bias+pose_bias+remainder, 'first_order_field': first,
            'half_second_order_field': second, 'joint_field': joint, 'minimum_field_bound': float(field_bound),
            'quadrature_kernel_error': error, 'stored_field_and_gram_bytes': storage,
            'quadrature_order': order, 'backend': backend}


def perturbed_geometry(k, q, pose, angle, shift):
    pose = np.asarray(pose)
    rotated = np.asarray(k)@Rotation.from_rotvec(angle*pose[:, :3]).as_matrix()
    translations = shift*np.einsum('nqa,na->nq', q, pose[:, 3:])
    return rotated, translations


def pose_cell_forward(k, q, ctf, coefficients, box, noise_std, pose, angle, shift):
    """Independent nonlinear projections of a constant-cell continuous density."""
    rotated, translations = perturbed_geometry(k, q, pose, angle, shift)
    n, nq = np.asarray(ctf).shape
    values = cell_forward(rotated, ctf, coefficients, box, noise_std).reshape(n, 2*nq)
    complex_values = (values[:, :nq]+1j*values[:, nq:])*np.exp(-2j*np.pi*translations)
    return np.concatenate([complex_values.real, complex_values.imag], axis=1).ravel()
