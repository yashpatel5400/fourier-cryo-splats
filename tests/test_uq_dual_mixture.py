import numpy as np
from fourier_splats.uq_dual_mixture import closest_support_mixture


def test_support_mixture_agrees_with_exact_segment_projection():
    rng = np.random.default_rng(609901)
    for _ in range(8):
        directions = rng.normal(size=(2, 11)); target = rng.normal(size=11)
        delta = directions[1]-directions[0]
        t = np.clip((target-directions[0])@delta/(delta@delta), 0., 1.)
        fit = closest_support_mixture(directions, target, initial=[.3, .7])
        np.testing.assert_allclose(fit['mixture'], [1-t, t], atol=2e-8)
        np.testing.assert_allclose(fit['distance'], np.linalg.norm(target-directions[0]-t*delta), rtol=1e-12)


def test_redundant_modes_preserve_feasibility_and_never_worsen_support():
    directions = np.array([[1., 0., 0.], [0., 1., 0.], [1., 0., 0.], [0., 0., 1.]])
    fit = closest_support_mixture(directions, np.full(3, 1/3), initial=[1., 0., 0., 0.])
    assert fit['distance'] < 1e-7
    assert fit['distance'] <= fit['initial_distance']
    assert np.all(fit['mixture'] >= 0)
    np.testing.assert_allclose(fit['mixture'].sum(), 1.)
