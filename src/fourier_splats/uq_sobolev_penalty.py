"""Convex derivative-tensor penalties for future pose-aware weight design.

No empirical optimization is performed here. Symmetric enclosing domains make
real and imaginary Fourier-coordinate Gram blocks orthogonal. Positive diagonal
pads preserve convexity in ordinary floating point; they are not interval bounds.
"""
import math
import numpy as np
from .uq_ball_remainder import particle_ball_remainder, cube_fourier_kernel, ball_fourier_kernel
from .uq_higher_remainder import partial_bell


def derivative_gram_blocks(k, transfer, radius, orders=(1, 2, 3, 4), domain='cube'):
    """Two real Gram blocks per derivative order, including transfer factors."""
    k = np.asarray(k, float); transfer = np.asarray(transfer, float)
    if k.ndim != 2 or k.shape != (len(transfer), 3) or not np.isfinite(k).all() or not np.isfinite(transfer).all():
        raise ValueError('Matching finite frequencies and transfer factors required')
    if not orders or any(not isinstance(m, (int, np.integer)) or m < 0 for m in orders):
        raise ValueError('Nonnegative derivative orders required')
    if domain not in ('cube', 'ball'):
        raise ValueError('Choose cube or ball')
    kernel = cube_fourier_kernel if domain == 'cube' else ball_fourier_kernel
    difference = kernel(k[:, None]-k[None, :], radius)
    addition = kernel(k[:, None]+k[None, :], radius)
    dot = k@k.T; transfer_product = transfer[:, None]*transfer[None, :]
    blocks = []; records = []
    for m in orders:
        common = .5*(2*np.pi)**(2*m)*dot**m*transfer_product
        pair = []
        for sign in (1, -1):
            gram = common*(difference+sign*(-1)**m*addition)
            gram = .5*(gram+gram.T)
            norm = float(np.linalg.norm(gram, ord=np.inf))
            minimum = float(np.linalg.eigvalsh(gram)[0])
            if minimum < -1e-10*max(norm, np.finfo(float).tiny):
                raise FloatingPointError('Substantially indefinite Fourier derivative Gram')
            # This is a declared heuristic arithmetic safeguard. In exact
            # arithmetic the integral Gram is PSD and needs no diagonal pad.
            pad = max(0., -minimum)+128*np.finfo(float).eps*len(k)*norm
            gram.flat[::len(k)+1] += pad
            pair.append(gram)
            records.append({'order': int(m), 'coordinate': 'real' if sign == 1 else 'imaginary',
                            'minimum_unpadded_eigenvalue': minimum, 'diagonal_pad': pad,
                            'infinity_norm': norm})
        blocks.append(pair)
    return np.asarray(blocks), records


class SobolevRemainderPenalty:
    """Convex degree-d field remainder and subgradient in real estimator weights.

    The real-arithmetic objective is a positive sum of Hilbert norms plus an
    amplitude residual penalty. It uses fixed geometry and joint pose balls.
    """
    def __init__(self, k, q, transfer, angle, shift, degree=3, domain='cube'):
        self.k = np.asarray(k, float); self.q = np.asarray(q, float)
        self.transfer = np.asarray(transfer, float)
        if self.transfer.ndim != 2:
            raise ValueError('Particle by frequency transfer required')
        self.n, self.nq = self.transfer.shape
        if self.k.shape != (self.n, self.nq, 3) or self.q.shape != (self.n, self.nq, 2):
            raise ValueError('Matching particle geometry required')
        if not isinstance(degree, (int, np.integer)) or not 2 <= degree <= 5:
            raise ValueError('Taylor degree between two and five required')
        self.order = degree+1
        angles = np.broadcast_to(angle, (self.n,)); shifts = np.broadcast_to(shift, (self.n,))
        self.blocks = np.empty((self.n, self.order, 2, self.nq, self.nq))
        self.multipliers = np.empty((self.n, self.order)); self.residual = np.empty((self.n, self.nq))
        self.diagnostics = []
        for i in range(self.n):
            a, s = angles[i], shifts[i]
            base = particle_ball_remainder(self.k[i], self.q[i], np.zeros(self.nq, complex), a, s, domain=domain)
            path = np.array([base['path_speed_bound']]+[a**j*np.sqrt(3)/2 for j in range(2, self.order+1)])
            self.multipliers[i] = partial_bell(path)[self.order, 1:]/math.factorial(self.order)
            self.blocks[i], diagnostics = derivative_gram_blocks(self.k[i], self.transfer[i], base['domain_radius'],
                orders=tuple(range(1, self.order+1)), domain=domain)
            embedding = np.linalg.lstsq(self.q[i], self.k[i], rcond=None)[0]
            eps = 2*np.pi*s*np.linalg.norm(self.q[i]-self.k[i]@np.linalg.pinv(embedding), axis=1)
            phase = path[:, None]*(2*np.pi*np.linalg.norm(self.k[i], axis=1))[None]
            complete = partial_bell(phase).sum(axis=1)
            residual = eps*complete[self.order]
            for j in range(1, self.order+1):
                residual += math.comb(self.order, j)*eps**j*complete[self.order-j]
            self.residual[i] = abs(self.transfer[i])*residual/math.factorial(self.order)
            self.diagnostics.append({'particle': i, 'domain_radius': base['domain_radius'], 'grams': diagnostics})

    def value_gradient(self, weights):
        w = np.asarray(weights, float).reshape(self.n, 2, self.nq)
        if not np.isfinite(w).all():
            raise ValueError('Finite weights required')
        actions = np.einsum('ndcqr,ncr->ndcq', self.blocks, w)
        squared = np.einsum('ndcq,ncq->nd', actions, w)
        if np.min(squared) < -1e-12:
            raise FloatingPointError('Negative PSD quadratic form')
        norms = np.sqrt(np.maximum(squared, 0.))
        factors = np.divide(self.multipliers, norms, out=np.zeros_like(norms), where=norms > 0)
        gradient = np.einsum('nd,ndcq->ncq', factors, actions)
        amplitude = np.hypot(w[:, 0], w[:, 1])
        ratios = np.divide(self.residual, amplitude, out=np.zeros_like(amplitude), where=amplitude > 0)
        gradient += ratios[:, None, :]*w
        value = float(np.sum(self.multipliers*norms)+np.sum(self.residual*amplitude))
        return value, gradient.ravel()
