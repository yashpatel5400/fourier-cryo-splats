"""Optional numerical preconditioner; does not change the inference operator."""
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse.linalg import LinearOperator


def diagonal_low_rank_inverse(factor, diagonal, ridge):
    """Inverse of FF' + diag(diag(G-FF')) + ridge I by Woodbury.

    A pivoted-Cholesky factor approximates G. Retaining its unrepresented
    diagonal avoids preconditioning every unresolved direction by ridge alone.
    This is a standard diagonal-plus-low-rank construction, not a new solver.
    """
    factor, diagonal = np.asarray(factor, float), np.asarray(diagonal, float)
    if (factor.ndim != 2 or diagonal.shape != (len(factor),) or ridge <= 0
            or not np.isfinite(ridge) or not np.isfinite(factor).all()
            or not np.isfinite(diagonal).all()):
        raise ValueError('Finite matching factor, diagonal and positive ridge required')
    residual = diagonal-np.sum(factor**2, axis=1)
    pad = 100*np.finfo(float).eps*max(1, factor.shape[1])*max(1., float(np.max(abs(diagonal))))
    if residual.min() < -pad:
        raise ValueError('Factor diagonal exceeds the original Gram diagonal')
    d = ridge+np.maximum(0., residual)
    weighted = factor/d[:, None]
    core = np.eye(factor.shape[1])+factor.T@weighted
    chol = cho_factor(core, lower=True, check_finite=False)
    def apply(x):
        return x/d-weighted@cho_solve(chol, factor.T@(x/d), check_finite=False)
    return LinearOperator((len(d), len(d)), matvec=apply, dtype=float)


def corrected_gram_preconditioner(gram, ridge):
    diagonal = gram.diagonal.reshape(gram.n, 2*gram.nq)
    blocks = [diagonal_low_rank_inverse(factor, d.ravel(), ridge)
        for (factor, _), d in zip(gram.preconditioner_factors,
            [diagonal[:, :gram.nq], diagonal[:, gram.nq:]])]
    if len(blocks) != 2:
        raise ValueError('Real and imaginary low-rank factors required')
    def apply(x):
        matrix = x.reshape(gram.n, 2*gram.nq)
        return gram.pack(blocks[0]@matrix[:, :gram.nq].ravel(),
                         blocks[1]@matrix[:, gram.nq:].ravel())
    return LinearOperator(gram.shape, matvec=apply, dtype=float)
