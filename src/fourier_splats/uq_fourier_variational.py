"""Hermitian Fourier-grid Gaussian VI for a fixed linear acquisition model.

Implements the Gaussian variational objective of Ullrich et al. (UAI/PMLR 115,
2020, Eq. 14) analytically for known poses/noise. CTFs, orthonormal real
coordinates, explicit prior scales and linear density functionals are project
extensions. No third-party implementation is copied. This is not joint pose
inference and does not reproduce the authors' complete experimental pipeline.
"""
import numpy as np
from scipy.sparse import coo_matrix
from .uq_continuous import gaussian_target_integrals, cell_forward


class HermitianTrilinearOperator:
    """Sparse trilinear interpolation of conjugate-symmetric Fourier values.

    An odd grid avoids ambiguous Nyquist coordinates. Parameters are the
    coefficients of 1, sqrt(2) cos(2*pi*f*x), sqrt(2) sin(2*pi*f*x), an
    orthonormal real Fourier basis on the unit cube. This fixes the prior and
    noise normalization. Interpolation is a model approximation, not the exact
    continuous transform at noninteger rotated frequencies.
    """
    def __init__(self, k, ctf, noise, box=33):
        k, ctf = np.asarray(k, float), np.asarray(ctf, float)
        if box < 3 or box % 2 != 1 or noise <= 0 or not np.isfinite(noise):
            raise ValueError('Positive noise and odd Fourier-grid size required')
        if k.ndim != 3 or k.shape[-1] != 3 or ctf.shape != k.shape[:2]:
            raise ValueError('Particle-frequency arrays required')
        self.box = box; self.n, self.nq = ctf.shape; h = box//2
        axis = np.arange(-h, h+1); z, y, x = np.meshgrid(axis, axis, axis, indexing='ij')
        all_f = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
        positive = self._sign(all_f) > 0
        self.frequencies = all_f[positive]
        lookup = np.full((box,)*3, -1, dtype=np.int32)
        f = self.frequencies
        lookup[f[:, 2]+h, f[:, 1]+h, f[:, 0]+h] = np.arange(len(f))
        points = k.reshape(-1, 3); lo = np.floor(points).astype(int); frac = points-lo
        if np.max(abs(points)) > h:
            raise ValueError('Observations outside the Fourier interpolation cube')
        transfer = ctf.ravel()/noise; count = len(points); index = np.arange(count)
        real_rows = (index//self.nq)*(2*self.nq)+index % self.nq
        rows, columns, values = [], [], []
        for mask in range(8):
            corner = np.array([(mask >> a) & 1 for a in range(3)])
            vertices = lo+corner; mass = np.prod(np.where(corner, frac, 1-frac), axis=1)
            active = np.flatnonzero(mass != 0.)
            v = vertices[active]; sign = self._sign(v); canonical = v*sign[:, None]
            ids = lookup[canonical[:, 2]+h, canonical[:, 1]+h, canonical[:, 0]+h]
            nonzero = sign != 0; coefficients = mass[active]*transfer[active]
            rows.append(real_rows[active]); columns.append(np.where(nonzero, 1+2*ids, 0))
            values.append(coefficients*np.where(nonzero, 1/np.sqrt(2), 1.))
            rows.append(real_rows[active[nonzero]]+self.nq)
            columns.append(2+2*ids[nonzero])
            values.append(-sign[nonzero]*coefficients[nonzero]/np.sqrt(2))
        self.matrix = coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
                                 shape=(2*count, box**3)).tocsr()
        self.matrix.sum_duplicates(); self.matrix.eliminate_zeros()
        self.gram_diagonal = np.asarray(self.matrix.power(2).sum(axis=0)).ravel()
        self.shape = self.matrix.shape

    @staticmethod
    def _sign(f):
        return np.where(f[:, 2] != 0, np.sign(f[:, 2]),
                        np.where(f[:, 1] != 0, np.sign(f[:, 1]), np.sign(f[:, 0]))).astype(int)

    def target(self, centers, signs, width):
        all_f = np.concatenate([np.zeros((1, 3)), self.frequencies])
        ft, full_norm2 = gaussian_target_integrals(all_f, centers, signs, width)
        ell = np.empty(self.shape[1]); ell[0] = ft[0].real
        ell[1::2] = np.sqrt(2)*ft[1:].real; ell[2::2] = -np.sqrt(2)*ft[1:].imag
        return ell, full_norm2

    def project_cells(self, rho, box):
        frequencies = np.concatenate([np.zeros((1, 3)), self.frequencies])[None]
        value = cell_forward(frequencies, np.ones(frequencies.shape[:2]), rho, box, 1.)
        count = frequencies.shape[1]
        ft = value[:count]+1j*value[count:]
        x = np.empty(self.shape[1]); x[0] = ft[0].real
        x[1::2] = np.sqrt(2)*ft[1:].real; x[2::2] = -np.sqrt(2)*ft[1:].imag
        return x


def diagonal_variational_variances(gram_diagonal, prior_sd):
    """Exact coordinate-variance optimum of the Gaussian negative ELBO.

    The mean equals the full Gaussian posterior mean, while q variance_j is
    1/precision_jj, not (precision^-1)_jj. A general target's credible variance
    is ell^T diag(1/precision_jj) ell. It is not a frequentist error bound.
    """
    diagonal = np.asarray(gram_diagonal, float)
    if prior_sd <= 0 or not np.isfinite(prior_sd) or (diagonal < 0).any():
        raise ValueError('Positive prior scale and nonnegative Gram diagonal required')
    return 1/(diagonal+1/prior_sd**2)
