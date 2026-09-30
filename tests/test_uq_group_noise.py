import numpy as np
import pytest
from fourier_splats.uq_group_noise import grouped_estimator_variance_upper, cell_density_distance
from fourier_splats.uq_noise_calibration import estimator_variance_upper


def test_singleton_groups_reproduce_independent_bound():
    rng = np.random.default_rng(371)
    y = rng.normal(size=(40, 7)); v = rng.normal(size=(5, 7))
    a = estimator_variance_upper(y, v)
    b = grouped_estimator_variance_upper(y, v, np.arange(5))
    assert b['noise_sd_upper'] == a['noise_sd_upper']


def test_perfect_dependence_attains_group_size_factor():
    rng = np.random.default_rng(372); y = rng.normal(size=(40, 7))
    v = np.tile(rng.normal(size=7), (5, 1))
    independent = estimator_variance_upper(y, v)
    grouped = grouped_estimator_variance_upper(y, v, np.zeros(5))
    np.testing.assert_allclose(grouped['estimator_variance_upper'], 5*independent['estimator_variance_upper'])
    # Every particle sharing the same vector noise attains the CS variance bound.
    np.testing.assert_allclose(np.linalg.norm(v.sum(axis=0))**2, 5*np.sum(v*v))


def test_arbitrary_joint_covariance_is_bounded():
    rng = np.random.default_rng(373); d = 4; n = 6
    common = rng.normal(size=(d, d)); sigma = common@common.T
    # Different orthogonal rotations of the same latent noise have equal
    # marginal covariance, with arbitrary cross-covariances across particles.
    factors = [common@np.linalg.qr(rng.normal(size=(d, d)))[0] for _ in range(n)]
    weights = rng.normal(size=(n, d))
    actual = np.linalg.norm(sum(w@a for w, a in zip(weights, factors)))**2
    envelope = n*sum(w@sigma@w for w in weights)
    assert actual <= envelope*(1+1e-12)


def test_nonnested_cells_match_exact_common_refinement():
    rng = np.random.default_rng(374); a = rng.normal(size=(3,)*3); b = rng.normal(size=(4,)*3)
    refine = lambda x, factor: x.repeat(factor, 0).repeat(factor, 1).repeat(factor, 2)/factor**1.5
    got, pad = cell_density_distance(a, 3, b, 4)
    np.testing.assert_allclose(got, np.linalg.norm(refine(a, 4)-refine(b, 3)), rtol=1e-13)
    assert pad > 0


def test_mismatched_groups_rejected():
    with pytest.raises(ValueError, match='label'):
        grouped_estimator_variance_upper(np.ones((5, 3)), np.ones((4, 3)), ['a'])
