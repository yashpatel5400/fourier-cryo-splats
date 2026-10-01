import numpy as np
from fourier_splats.uq_end_to_end import dephase_observations, observed_interval, binomial_interval


def test_dephase_is_orthogonal_and_removes_known_translation():
    rng = np.random.default_rng(33)
    q = rng.normal(size=(4, 7, 2)); shift = rng.normal(size=(4, 2)); y = rng.normal(size=(4, 14))
    shifted = dephase_observations(y, q, -shift, 200).reshape(4, 14)
    np.testing.assert_allclose(dephase_observations(shifted, q, shift, 200), y.ravel(), atol=1e-14)
    np.testing.assert_allclose(np.sum(shifted**2, axis=1), np.sum(y**2, axis=1), atol=1e-13)


def test_fallback_does_not_hide_raw_failure_and_retains_center():
    r = observed_interval(50., 3., 1., 2., 0., True)
    assert not r['raw_covered'] and r['covered'] and r['fallback']
    assert r['center'] == 1. and r['raw_center'] == 50.
    assert not r['correct_sign_exclusion']  # Zero truth is not a successful sign conclusion.
    raw = observed_interval(50., 3., 1., 2., 0., False)
    assert not raw['covered'] and not raw['fallback']


def test_exact_binomial_interval_includes_boundary_and_expected_range():
    lo, hi = binomial_interval(190, 200)
    assert .90 < lo < .92 and .97 < hi < .98
    assert binomial_interval(0, 200)[0] == 0.
    assert binomial_interval(200, 200)[1] == 1.
