"""Local Gaussian nuisance baseline; no uniform nonlinear pose guarantee.

The Jacobian differentiates the same trilinear Fourier model as the baseline,
at its fixed pilot. Density/pose cross terms and second derivatives are omitted.
The Gaussian pose prior is supplied, not estimated or calibrated here.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, aslinearoperator


def trilinear_pose_jacobian(operator, k, q, ctf, noise, coefficients):
    """Per-particle derivatives for row rotations and detector translations.

    Columns are three radians and two field-fraction translations. The Fourier
    interpolation derivative exists away from its integer knots; reject knots
    rather than silently choosing a one-sided derivative for a claimed Jacobian.
    """
    k, q, ctf = map(lambda x: np.asarray(x, float), (k, q, ctf))
    x = np.asarray(coefficients, float); n, nq = ctf.shape; h = operator.box//2
    if k.shape != (n, nq, 3) or q.shape != (n, nq, 2) or x.shape != (operator.shape[1],):
        raise ValueError('Matching finite-model coefficients and particle geometry required')
    if noise <= 0 or not np.isfinite(noise) or not all(np.isfinite(a).all() for a in [k, q, ctf, x]):
        raise ValueError('Positive noise and finite inputs required')
    distance = float(np.min(abs(k-np.rint(k))))
    if distance <= 1e-10 or np.max(abs(k)) >= h:
        raise ValueError('Pose derivative requires interior non-knot interpolation points')
    grid = np.zeros((operator.box,)*3, complex); grid[h, h, h] = x[0]
    f = operator.frequencies; values = (x[1::2]-1j*x[2::2])/np.sqrt(2)
    grid[tuple((f+h).T[::-1])] = values
    grid[tuple((-f+h).T[::-1])] = values.conj()
    lo = np.floor(k).astype(int); frac = k-lo
    value = np.zeros((n, nq), complex); gradient = np.zeros((n, nq, 3), complex)
    for mask in range(8):
        corner = np.array([(mask >> a) & 1 for a in range(3)])
        vertices = lo+corner+h
        selected = grid[vertices[..., 2], vertices[..., 1], vertices[..., 0]]
        factors = np.where(corner, frac, 1-frac)
        value += selected*np.prod(factors, axis=-1)
        for axis in range(3):
            gradient[..., axis] += selected*(2*corner[axis]-1)*np.prod(np.delete(factors, axis, axis=-1), axis=-1)
    derivative = np.empty((n, nq, 5), complex)
    for axis in range(3):
        # k [e_axis]_cross equals k cross e_axis for row-vector rotations.
        tangent = np.cross(k, np.eye(3)[axis])
        derivative[..., axis] = np.sum(gradient*tangent, axis=-1)
    derivative[..., 3:] = -2j*np.pi*q*value[..., None]
    derivative *= (ctf/noise)[..., None]
    return np.concatenate([derivative.real, derivative.imag], axis=1), {
        'minimum_distance_to_integer_knot': distance,
        'scope': 'Derivative of trilinear pilot signal, not continuous-density or nonlinear joint inference.'}


class GaussianNuisanceWhitening:
    """Apply (I+J J^T)^(-1/2) with independent small particle blocks."""
    def __init__(self, scaled_jacobian):
        j = np.asarray(scaled_jacobian, float)
        if j.ndim != 3 or not np.isfinite(j).all():
            raise ValueError('Finite particle by measurement by nuisance Jacobian required')
        self.j = j; self.n, self.m, self.d = j.shape
        gram = j.swapaxes(-1, -2)@j
        eigenvalues, vectors = np.linalg.eigh(gram)
        if np.min(eigenvalues) < -1e-12*max(1., np.max(eigenvalues)):
            raise FloatingPointError('Negative nuisance Gram eigenvalue')
        values = np.maximum(eigenvalues, 0.)
        self.factors = j@vectors
        root = np.sqrt(1+values)
        self.multipliers = -1/(root*(root+1))

    def apply(self, vector):
        v = np.asarray(vector, float).reshape(self.n, self.m)
        small = np.einsum('nmd,nm->nd', self.factors, v)
        return (v+np.einsum('nmd,nd->nm', self.factors, self.multipliers*small)).ravel()

    def covariance_apply(self, vector):
        v = np.asarray(vector, float).reshape(self.n, self.m)
        small = np.einsum('nmd,nm->nd', self.j, v)
        return (v+np.einsum('nmd,nd->nm', self.j, small)).ravel()

    def whiten_operator(self, operator):
        a = aslinearoperator(operator)
        if a.shape[0] != self.n*self.m:
            raise ValueError('Nuisance blocks and acquisition rows disagree')
        return LinearOperator(a.shape, matvec=lambda v: self.apply(a@v),
            rmatvec=lambda v: a.rmatvec(self.apply(v)), dtype=float)
