"""Cubic pose lift with twenty continuous spatial moments.

Positive block scaling and quadrature integration are separate operations.
This module does not calibrate experimental poses or density-class radii.
"""
import itertools
import math
import numpy as np
import finufft
from scipy.special import spherical_jn, ndtr
from .uq_pose_operator import PolynomialPoseFieldOperator
from .uq_joint_bias import MONOMIALS as QUADRATIC_MONOMIALS
from .uq_continuous_pose import cube_quadrature, polynomial_kernel_error
from .uq_continuous_quadrature import QuadratureObservationGram

SPATIAL = tuple(QUADRATIC_MONOMIALS)+tuple(
    tuple(t.count(i) for i in range(3))
    for t in itertools.combinations_with_replacement(range(3), 3))
LOOKUP = {p: j for j, p in enumerate(SPATIAL)}
POSE = tuple(tuple(t.count(i) for i in range(5))
             for degree in (1, 2, 3)
             for t in itertools.combinations_with_replacement(range(5), degree))
BLOCK_SIZES = (5, 15, 35)
NORMALIZATION = np.array([math.sqrt(math.factorial(sum(p))/math.prod(math.factorial(v) for v in p)) for p in POSE])


def multiply_polynomials(left, right):
    """Multiply spatial polynomials whose total product degree is at most 3."""
    out = np.zeros(np.broadcast_shapes(left.shape[:-1], right.shape[:-1])+(20,), dtype=np.result_type(left, right))
    for a, pa in enumerate(SPATIAL[:left.shape[-1]]):
        for b, pb in enumerate(SPATIAL[:right.shape[-1]]):
            total = tuple(x+y for x, y in zip(pa, pb))
            if sum(total) <= 3:
                out[..., LOOKUP[total]] += left[..., a]*right[..., b]
    return out


def pose_lift(u):
    u = np.asarray(u, float)
    if u.shape[-1] != 5:
        raise ValueError('Five pose coordinates required')
    return np.stack([scale*np.prod(u**np.array(alpha), axis=-1) for alpha, scale in zip(POSE, NORMALIZATION)], axis=-1)


class CubicPoseFieldOperator(PolynomialPoseFieldOperator):
    """Degree-major columns, with particle-major columns inside each degree.

    Geometry coefficients are cached once. This bounded diagnostic is intended
    for small cohorts; storage of this explicit coefficient tensor is reported.
    It does not store the much larger quadrature-by-particle field matrix.
    """
    def __init__(self, *args, **kwargs):
        if kwargs.get('column_denominators') is not None:
            raise ValueError('Set cubic scales explicitly after initialization')
        super().__init__(*args, **kwargs)
        self.shape = (len(self.xyz), 55*self.n)
        self.denominators = np.ones(55*self.n)
        self.monomials = np.stack([np.prod(self.xyz**np.array(beta), axis=1) for beta in SPATIAL])
        self.polynomials = np.empty((self.n, self.nq, 55, 20), complex)
        for index in self._blocks():
            self.polynomials[index] = self._cubic_coefficients(index)

    def _cubic_coefficients(self, index):
        b = self._phase_coefficients(index)
        out = np.zeros((*b.shape[:2], 55, 20), complex)
        out[..., :5, :4] = 1j*b
        phi2 = {}
        for j, alpha in enumerate(POSE[5:20]):
            axes = [a for a, count in enumerate(alpha) for _ in range(count)]
            a, bb = axes; factor = .5 if a == bb else 1.
            v = np.zeros((*b.shape[:2], 4))
            v[..., 1:] = factor*self._curvature(index, a, bb)
            phi2[alpha] = v
            out[..., 5+j, :4] = 1j*v
            out[..., 5+j, :] -= factor*multiply_polynomials(b[..., a, :], b[..., bb, :])
        for j, alpha in enumerate(POSE[20:]):
            col = out[..., 20+j, :]
            axes = [a for a, count in enumerate(alpha) for _ in range(count)]
            product = multiply_polynomials(b[..., axes[0], :], b[..., axes[1], :])[..., :10]
            col -= 1j*multiply_polynomials(product, b[..., axes[2], :])/math.prod(math.factorial(v) for v in alpha)
            for a in range(5):
                if alpha[a]:
                    other = list(alpha); other[a] -= 1
                    col -= multiply_polynomials(b[..., a, :], phi2[tuple(other)])
            for a in range(3):
                for c in range(3):
                    powers = np.zeros(5, int); powers[a] += 1; powers[c] += 2
                    if tuple(powers) == alpha:
                        col[..., :4] -= 1j*self.angle[index, None, None]**2*b[..., a, :]/6
        return out/NORMALIZATION[None, None, :, None]

    def _to_particle_columns(self, vector):
        vector = np.asarray(vector, float).reshape(55*self.n)
        ends = np.cumsum([0]+[size*self.n for size in BLOCK_SIZES])
        return np.concatenate([vector[a:b].reshape(self.n, size)
                               for a, b, size in zip(ends[:-1], ends[1:], BLOCK_SIZES)], axis=1)

    @staticmethod
    def _from_particle_columns(value):
        return np.concatenate([value[:, :5].ravel(), value[:, 5:20].ravel(), value[:, 20:].ravel()])

    def _transform_forward(self, strengths):
        if self.backend == 'nufft':
            return finufft.nufft3d3(*self.source, np.ascontiguousarray(strengths), *self.destination,
                                    isign=1, eps=self.eps, nthreads=self.nthreads)
        output = np.empty((20, len(self.xyz)), complex)
        for start in range(0, len(self.xyz), 1024):
            phase = np.exp(2j*np.pi*self.frequency_points@self.xyz[start:start+1024].T)
            output[:, start:start+phase.shape[1]] = strengths@phase
        return output

    def _transform_adjoint(self, strengths):
        if self.backend == 'nufft':
            return finufft.nufft3d3(*self.destination, np.ascontiguousarray(strengths), *self.source,
                                    isign=1, eps=self.eps, nthreads=self.nthreads)
        output = np.empty((20, len(self.frequency_points)), complex)
        for start in range(0, len(self.frequency_points), 1024):
            phase = np.exp(2j*np.pi*self.xyz@self.frequency_points[start:start+1024].T)
            output[:, start:start+phase.shape[1]] = strengths@phase
        return output

    def matvec(self, vector):
        columns = self._to_particle_columns(np.asarray(vector)/self.denominators)
        polynomial = np.einsum('nqap,na->nqp', self.polynomials, columns)
        strengths = self.c[..., None]*polynomial
        field = self._transform_forward(strengths.reshape(-1, 20).T)
        return self.sqrt_quad*np.sum(self.monomials*field.real, axis=0)

    def pair_moments(self, moments):
        moments = np.asarray(moments)
        if moments.shape != (*self.k.shape[:-1], 20):
            raise ValueError('Twenty moments per frequency required')
        value = np.einsum('nqap,nqp->na', self.polynomials, moments*self.c[..., None]).real
        return self._from_particle_columns(value)/self.denominators

    def rmatvec(self, vector):
        vector = np.asarray(vector, float).reshape(len(self.xyz))
        transformed = self._transform_adjoint((self.monomials*(self.sqrt_quad*vector)[None]).astype(complex))
        return self.pair_moments(transformed.T.reshape(self.n, self.nq, 20))

    def weight_gradient(self, *args, **kwargs):
        raise NotImplementedError('This fixed-weight cubic diagnostic does not optimize weights')

    def coefficient_bound_matrix(self):
        den = self._to_particle_columns(self.denominators)
        return np.sum(abs(self.polynomials), axis=-1).transpose(0, 2, 1)*abs(self.transfer)[:, None, :]/den[:, :, None]

    def establish_quadrature_design_scaling(self, order=12, callback=None):
        if not isinstance(order, (int, np.integer)) or order < 2:
            raise ValueError('Integer design quadrature order >=2 required')
        xyz, qw = cube_quadrature(order)
        monomials = np.stack([np.prod(xyz**np.array(beta), axis=1) for beta in SPATIAL])
        squared = np.empty((self.n, 3))
        for i in range(self.n):
            phase = np.exp(2j*np.pi*self.k[i]@xyz.T)*self.c[i, :, None]
            fields = np.einsum('qap,qr,pr->ra', self.polynomials[i], phase, monomials, optimize=True).real
            for j, (a, b) in enumerate(((0, 5), (5, 20), (20, 55))):
                squared[i, j] = np.sum(qw[:, None]*fields[:, a:b]**2)
            if callback and ((i+1) % 32 == 0 or i+1 == self.n):
                callback({'cubic_design_particles': i+1, 'total': self.n})
        safe = np.where(squared > 0, np.sqrt(squared), 1.)
        if not np.isfinite(safe).all():
            raise FloatingPointError('Nonfinite cubic scales')
        scales = safe.T.ravel()
        self.denominators = np.sqrt(np.concatenate([np.repeat(safe[:, j], size) for j, size in enumerate(BLOCK_SIZES)]))
        return {'group_scales': scales, 'sum_group_scales': float(scales.sum()), 'design_order': order,
                'scaling_method': 'Positive sampled block norms used only to choose scales; not integration certificates'}

    def establish_group_scaling(self, *args, **kwargs):
        raise NotImplementedError('Use explicit cubic design scaling and a separate padded spectral audit')

    def establish_coefficient_scaling(self, *args, **kwargs):
        raise NotImplementedError('Use the declared cubic design scaling')


def normalized_cell_moments_cubic(k, box):
    """Analytic moments via Legendre expansions and spherical Bessel functions."""
    x = np.pi*np.asarray(k, float)/box; a = 1/(2*box)
    j = [spherical_jn(n, abs(x))*np.where(x < 0, (-1)**n, 1) for n in range(4)]
    return np.stack([j[0], 1j*a*j[1], a*a*(j[0]-2*j[2])/3,
                     1j*a**3*(3*j[1]-2*j[3])/5], axis=-1)


def cell_moments_cubic(k, coefficients, box, nthreads=1):
    k = np.asarray(k, float); coef = np.asarray(coefficients, float)
    if box < 2 or box % 2 or coef.size != box**3 or k.shape[-1] != 3:
        raise ValueError('Matching even constant-cell grid and 3D frequencies required')
    if not np.isfinite(k).all() or not np.isfinite(coef).all():
        raise ValueError('Finite cell-moment inputs required')
    coordinate = (np.arange(box)+.5)/box-.5
    z, y, x = np.meshgrid(coordinate, coordinate, coordinate, indexing='ij'); xyz = [x, y, z]
    strengths = np.stack([coef.reshape((box,)*3)*np.prod(np.stack([xyz[a]**p for a, p in enumerate(beta)]), axis=0) for beta in SPATIAL]).astype(complex)
    points = [np.ascontiguousarray(2*np.pi*k.reshape(-1, 3)[:, axis]/box) for axis in (2, 1, 0)]
    transform = finufft.nufft3d2(*points, strengths, isign=1, eps=1e-12, nthreads=nthreads)
    transform = transform.T.reshape((*k.shape[:-1], 20))*np.exp(1j*np.pi*k.sum(axis=-1)/box)[..., None]/box**1.5
    local = normalized_cell_moments_cubic(k, box); result = np.zeros_like(transform)
    for j, beta in enumerate(SPATIAL):
        for gamma in itertools.product(*(range(p+1) for p in beta)):
            factor = np.ones(k.shape[:-1], complex)
            for axis, (b, g) in enumerate(zip(beta, gamma)):
                factor *= math.comb(b, g)*local[..., axis, b-g]
            result[..., j] += transform[..., LOOKUP[gamma]]*factor
    return result


def gaussian_moments_cubic(k, centers, signs, width):
    k = np.asarray(k, float); centers = np.asarray(centers, float); sigma = float(width)
    if sigma <= 0 or centers.shape != (len(signs), 3) or k.shape[-1] != 3:
        raise ValueError('Valid Gaussian target required')
    omega = 2*np.pi*k; result = np.zeros((*k.shape[:-1], 20), complex)
    for mu, sign in zip(centers, signs):
        shift = mu+1j*sigma*sigma*omega
        zero = np.exp(1j*omega*mu-.5*sigma*sigma*omega**2)*(ndtr((.5-mu)/sigma-1j*sigma*omega)-ndtr((-.5-mu)/sigma-1j*sigma*omega))
        upper = np.exp(-(.5-mu)**2/(2*sigma*sigma)+.5j*omega)/(np.sqrt(2*np.pi)*sigma)
        lower = np.exp(-(-.5-mu)**2/(2*sigma*sigma)-.5j*omega)/(np.sqrt(2*np.pi)*sigma)
        moments = [zero]
        for degree in range(1, 4):
            value = shift*moments[-1]-sigma*sigma*(.5**(degree-1)*upper-(-.5)**(degree-1)*lower)
            if degree > 1:
                value += (degree-1)*sigma*sigma*moments[-2]
            moments.append(value)
        for j, beta in enumerate(SPATIAL):
            result[..., j] += sign*np.prod(np.stack([moments[p][..., axis] for axis, p in enumerate(beta)]), axis=0)
    if not np.isfinite(result).all():
        raise FloatingPointError('Nonfinite Gaussian moments')
    return result


def cubic_residual_cross(op, weights, noise, centers, signs, width):
    gram = QuadratureObservationGram(op.k, op.ctf, noise, order=op.order, preconditioner_rank=0)
    gram.nthreads = op.nthreads
    target = op.pair_moments(gaussian_moments_cubic(op.k, centers, signs, width))
    cross = target-op.rmatvec(op.sqrt_quad*gram.field(weights))
    w = np.asarray(weights).reshape(op.n, 2*op.nq); amp = np.hypot(w[:, :op.nq], w[:, op.nq:])
    masses = np.einsum('naj,nj->na', op.coefficient_bound_matrix(), amp)
    pad = polynomial_kernel_error(op.k, op.order, degree=6)*np.sum(abs(gram.coefficients(weights)))*np.linalg.norm(masses)
    return {'norm_upper': float(np.linalg.norm(cross)+pad), 'norm_unpadded': float(np.linalg.norm(cross)),
            'quadrature_norm_pad': float(pad), 'cross_vector': cross}


def direct_cell_moments_cubic(k, coefficients, box, chunk=16):
    """Selected-frequency check using physical cells and 24-node local quadrature.

    The local quadrature is an independent numerical check on Bessel formulas,
    not a general integration certificate. No NUFFT is used here.
    """
    from numpy.polynomial.legendre import leggauss
    k = np.asarray(k, float); coef = np.asarray(coefficients, float).reshape(-1)
    if coef.size != box**3:
        raise ValueError('Matching cell coefficients required')
    centers = (np.indices((box,)*3).reshape(3, -1).T[:, ::-1]+.5)/box-.5
    nodes, weights = leggauss(24); offsets = nodes/(2*box)
    frequencies = k.reshape(-1, 3); result = np.empty((len(frequencies), 20), complex)
    for start in range(0, len(frequencies), chunk):
        frequency = frequencies[start:start+chunk]
        phase = np.exp(2j*np.pi*frequency@centers.T)
        local = np.stack([np.exp(2j*np.pi*frequency[..., None]*offsets)@(weights*offsets**degree/2) for degree in range(4)], axis=-1)
        translated = []
        for axis in range(3):
            translated.append([sum(math.comb(degree, j)*centers[None, :, axis]**(degree-j)*local[:, axis, j, None]
                                   for j in range(degree+1)) for degree in range(4)])
        for j, beta in enumerate(SPATIAL):
            moment = phase.copy()
            for axis, power in enumerate(beta):
                moment *= translated[axis][power]
            result[start:start+len(frequency), j] = moment@coef/box**1.5
    return result.reshape((*k.shape[:-1], 20))
