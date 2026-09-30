import numpy as np
from scipy.interpolate import RegularGridInterpolator
from fourier_splats.uq_fourier_variational import HermitianTrilinearOperator, diagonal_variational_variances
from fourier_splats.uq_baselines import gaussian_reference_operator


def test_hermitian_sparse_slice_matches_independent_complex_interpolation():
    rng = np.random.default_rng(961); k = rng.uniform(-2.8, 2.8, size=(4, 7, 3))
    ctf = rng.normal(size=(4, 7)); op = HermitianTrilinearOperator(k, ctf, .7, box=7)
    x = rng.normal(size=op.shape[1]); grid = np.zeros((7,)*3, complex); h = 3
    grid[h, h, h] = x[0]
    for f, a, b in zip(op.frequencies, x[1::2], x[2::2]):
        grid[tuple((f+h)[::-1])] = (a-1j*b)/np.sqrt(2)
        grid[tuple((-f+h)[::-1])] = (a+1j*b)/np.sqrt(2)
    rgi = RegularGridInterpolator((np.arange(-h, h+1),)*3, grid)
    pred = rgi(k[..., ::-1].reshape(-1, 3)).reshape(ctf.shape)*ctf/.7
    expected = np.concatenate([pred.real, pred.imag], axis=1).ravel()
    np.testing.assert_allclose(op.matrix@x, expected, atol=1e-14)
    np.testing.assert_allclose(np.sum(abs(grid)**2), x@x, rtol=1e-14)
    negative = HermitianTrilinearOperator(-k, ctf, .7, box=7).matrix@x
    np.testing.assert_allclose(negative.reshape(4, 14)[:, :7], pred.real, atol=1e-14)
    np.testing.assert_allclose(negative.reshape(4, 14)[:, 7:], -pred.imag, atol=1e-14)


def test_diagonal_elbo_stationarity_and_full_posterior_target_solution():
    rng = np.random.default_rng(962); k = rng.uniform(-.8, .8, (3, 4, 3))
    op = HermitianTrilinearOperator(k, np.ones((3, 4)), .3, box=3)
    a = op.matrix.toarray(); tau = .4; ell, _ = op.target([[.03, .02, 0]], [1], .1)
    precision = a.T@a+np.eye(a.shape[1])/tau**2
    variance = diagonal_variational_variances(op.gram_diagonal, tau)
    np.testing.assert_allclose(np.diag(precision)-1/variance, 0., atol=2e-14)
    y = rng.normal(size=a.shape[0]); mean = np.linalg.solve(precision, a.T@y)
    np.testing.assert_allclose(a.T@(a@mean-y)+mean/tau**2, 0., atol=2e-14)
    fit = gaussian_reference_operator(op.matrix, ell, tau, gram_diagonal=op.gram_diagonal)
    np.testing.assert_allclose(fit['weights']@y, ell@mean, rtol=1e-9)
    residual = a.T@fit['weights']-ell
    np.testing.assert_allclose(fit['weights']@fit['weights']+tau**2*(residual@residual),
                              ell@np.linalg.solve(precision, ell), rtol=1e-9)
