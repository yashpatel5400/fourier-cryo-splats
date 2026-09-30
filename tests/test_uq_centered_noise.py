import numpy as np
from fourier_splats.uq_centered_noise import centered_grouped_variance_upper
from fourier_splats.uq_noise_calibration import common_covariance_trace_upper


def test_orthogonal_contrasts_and_common_shift():
    rng = np.random.default_rng(951230)
    y = rng.normal(size=(18, 7)); v = rng.normal(size=(9, 7)); groups = np.arange(9)//2
    result, q, contrasted = centered_grouped_variance_upper(y, v, groups)
    np.testing.assert_allclose(q@q.T, np.eye(17), atol=5e-16)
    np.testing.assert_allclose(q@np.ones(18), 0, atol=5e-16)
    shifted, _, cy = centered_grouped_variance_upper(y+rng.normal(size=7)*30, v, groups)
    np.testing.assert_allclose(contrasted, cy, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(result['noise_sd_upper'], shifted['noise_sd_upper'], rtol=1e-12)
    sizes = np.bincount(groups)[groups]; weighted = v*np.sqrt(sizes[:, None])
    explicit_energy = np.sum(((y-y.mean(axis=0))@weighted.T)**2)
    np.testing.assert_allclose(result['calibration']['observed_total_energy'], explicit_energy, rtol=1e-13)
    assert result['calibration']['observations'] == 17


def test_unequal_mean_gaussian_covariance_and_bound():
    # Independent known-covariance validation, including arbitrary unequal means.
    rng = np.random.default_rng(951231); draws = 8000; n = 12; dim = 3
    factor = np.array([[1., 0., 0.], [.5, .3, 0.], [-.2, .1, .1]])
    sigma = factor@factor.T
    _, q, _ = centered_grouped_variance_upper(np.zeros((n, dim)), np.eye(dim), np.arange(dim), .05)
    means = rng.normal(size=(n, dim))*.1+np.array([5., -2., 1.])
    z = rng.normal(size=(draws, n, dim))@factor.T+means
    contrasts = np.einsum('ij,bjk->bik', q, z)
    centered_noise = contrasts-q@means
    np.testing.assert_allclose(np.einsum('bik,bil->kl', centered_noise, centered_noise)/(draws*(n-1)), sigma, atol=.017)
    # Central rank-one gives the loosest transform inequality; this is an
    # implementation check, not empirical evidence for archive noise models.
    fraction = common_covariance_trace_upper(np.zeros((n-1, dim)), .05)['lower_tail_fraction']
    upper = np.sum(contrasts**2, axis=(1, 2))/((n-1)*fraction)
    failure = np.mean(upper < np.trace(sigma))
    assert failure < .05 + 4*np.sqrt(.05*.95/draws)


def test_invalid_one_row_rejected():
    import pytest
    with pytest.raises(ValueError, match='At least two'):
        centered_grouped_variance_upper(np.zeros((1, 2)), np.ones((1, 2)), np.array([0]))
