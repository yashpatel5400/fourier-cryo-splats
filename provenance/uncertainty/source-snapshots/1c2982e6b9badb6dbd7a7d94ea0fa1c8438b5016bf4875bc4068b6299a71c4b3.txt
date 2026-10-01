import numpy as np
import pytest
from fourier_splats.uq_movie_diagnostics import radial_power, correlation, block_average


def test_real_fourier_shells_preserve_variance():
    rng = np.random.default_rng(55)
    x = rng.normal(size=(32, 32))+8
    d = radial_power(x)
    assert np.sum(d['shell_energy']) == pytest.approx(np.var(x), rel=1e-12)
    assert np.sum(d['modes']) == x.size
    assert d['shell_energy'][0] < 1e-25


def test_disjoint_and_shared_frame_differences_have_distinct_nulls():
    rng = np.random.default_rng(61)
    frames = rng.normal(size=(4, 512, 512))
    assert abs(correlation(frames[1]-frames[0], frames[3]-frames[2])) < .01
    assert correlation(frames[1]-frames[0], frames[2]-frames[1]) == pytest.approx(-.5, abs=.01)


def test_block_average_preserves_mean_without_interpolation():
    x = np.arange(64).reshape(8, 8)
    y = block_average(x, 2)
    assert y.shape == (4, 4)
    assert y.mean() == x.mean()
    assert y[0, 0] == 4.5
