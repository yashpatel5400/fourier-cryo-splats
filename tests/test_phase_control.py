import numpy as np
import pytest
from fourier_splats.uq_phase_control import cross_power_amplitude_interval


@pytest.mark.parametrize('sigma,n,mean', [(2., 16, .4), (.25, 4096, 1.01), (3., 1, -.5)])
def test_inversion_agrees_with_original_confidence_inequality(sigma, n, mean):
    interval = cross_power_amplitude_interval(mean, n, sigma)
    assert not interval['empty']
    x = np.log(40.)
    for s in np.linspace(0., 2*interval['upper']+1, 1001):
        accepted = abs(mean-s*s) <= 2*np.sqrt((sigma**4+s*s*sigma*sigma)*x/n)+sigma*sigma*x/n
        assert accepted == (interval['lower'] <= s <= interval['upper'])


def test_empty_set_is_not_silently_replaced_and_noise_is_checked():
    interval = cross_power_amplitude_interval(-100., 100, 1.)
    assert interval == dict(empty=True, lower=None, upper=None)
    with pytest.raises(ValueError):
        cross_power_amplitude_interval(1., 100, 0.)


def test_product_normal_moments_against_independent_gaussian_samples():
    rng = np.random.default_rng(781)
    s, sigma, count = 1.3, .7, 120000
    a = s+sigma*(rng.normal(size=count)+1j*rng.normal(size=count))
    b = s+sigma*(rng.normal(size=count)+1j*rng.normal(size=count))
    product = (a.conj()*b).real
    expected_variance = 2*(sigma**4+s*s*sigma*sigma)
    assert abs(product.mean()-s*s) < 5*np.sqrt(expected_variance/count)
    np.testing.assert_allclose(product.var(), expected_variance, rtol=.025)
    for t in [-.08, .08]:
        expected = np.exp(-np.log1p(-sigma**4*t*t)+s*s*sigma*sigma*t*t/(1-sigma*sigma*t))
        np.testing.assert_allclose(np.exp(t*(product-s*s)).mean(), expected, rtol=.002)
