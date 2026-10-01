import numpy as np
import pytest
from fourier_splats.uq_gaussian_preconditioner import diagonal_low_rank_inverse


@pytest.mark.parametrize('ridge', [.01, 1., 100.])
def test_inverse_matches_dense_diagonal_plus_low_rank(ridge):
    rng = np.random.default_rng(81)
    factor = rng.normal(size=(30, 7))
    remainder = rng.uniform(.02, 4, 30)
    diagonal = np.sum(factor**2, axis=1)+remainder
    matrix = factor@factor.T+np.diag(remainder+ridge)
    pre = diagonal_low_rank_inverse(factor, diagonal, ridge)
    b = rng.normal(size=30)
    np.testing.assert_allclose(pre@b, np.linalg.solve(matrix, b), rtol=1e-11, atol=1e-11)
    assert b@(pre@b) > 0
    with pytest.raises(ValueError):
        diagonal_low_rank_inverse(factor, diagonal-100, ridge)
