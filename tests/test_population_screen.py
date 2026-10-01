import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from scipy.stats import norm
from fourier_splats.uq_population_screen import (
    amplitude_tangent_information, gaussian_mixture_information,
    pseudo_true_fraction, histogram_law_weights, stratified_linear_summary)


def test_gaussian_information_against_direct_density_quadrature():
    for s in [.0001, .2, 1., 4.]:
        for t in [.2, .75]:
            d = np.sqrt(s)
            def integrand(y):
                a, b = norm.pdf(y-d), norm.pdf(y)
                return (a-b)**2 / (t*a+(1-t)*b)
            direct = quad(integrand, -12, 12+d, epsabs=1e-11)[0]
            actual = gaussian_mixture_information(np.array([s]), t, 128)[0]
            np.testing.assert_allclose(actual, direct, atol=2e-11, rtol=1e-9)


def test_population_root_against_expected_log_density_optimization():
    energy = np.array([.2, 2.])
    wa, wb = np.array([.1, .9]), np.array([.8, .2])
    p = .75
    def objective(t):
        total = 0.
        for s, a, b in zip(energy, wa, wb):
            d = np.sqrt(s)
            def integrand(y):
                la, lb = norm.logpdf(y-d), norm.logpdf(y)
                truth = p*a*np.exp(la)+(1-p)*b*np.exp(lb)
                return truth * logsumexp([np.log(t)+la, np.log1p(-t)+lb])
            total += quad(integrand, -12, 12+d, epsabs=1e-10)[0]
        return -total
    opt = minimize_scalar(objective, bounds=(.001, .999), method='bounded', options={'xatol':1e-11})
    result = pseudo_true_fraction(energy, wa, wb)
    np.testing.assert_allclose(result['fraction'], opt.x, atol=2e-8, rtol=0)
    assert result['quadrature_check_passed']
    # A shared view law cancels the misspecification effect, despite unequal s.
    shared = pseudo_true_fraction(energy, wa, wa)
    np.testing.assert_allclose(shared['fraction'], p, atol=1e-10)


def test_real_amplitude_tangent_by_explicit_least_squares():
    rng = np.random.default_rng(811)
    m = rng.normal(size=(8,7)) + 1j*rng.normal(size=(8,7))
    g = rng.normal(size=(8,7)) + 1j*rng.normal(size=(8,7))
    s, a, b = amplitude_tangent_information(m, g)
    for point, actual in [(m,a),(m-g,b)]:
        expected = []
        for mm, gg in zip(point,g):
            x, y = np.r_[mm.real,mm.imag], np.r_[gg.real,gg.imag]
            residual = y-x*np.linalg.lstsq(x[:,None],y,rcond=None)[0][0]
            expected.append(residual@residual)
        np.testing.assert_allclose(actual, expected, atol=1e-13)
    # An imaginary multiple is real-amplitude orthogonal, not complex-collinear.
    _, residual, _ = amplitude_tangent_information(np.ones((1,2),complex), 1j*np.ones((1,2)))
    assert residual[0] == 2


def test_stratified_weighting_and_conditional_sampling_error():
    bins = np.array([0,0,1,1,1])
    values = np.array([1.,3.,4.,5.,9.])
    counts = np.array([1,3])
    w = histogram_law_weights(bins,counts)
    np.testing.assert_allclose(w.sum(),1)
    expected = .25*2 + .75*6
    result = stratified_linear_summary(values,bins,counts)
    np.testing.assert_allclose(result['mean'],expected)
    variance = .25**2*np.var(values[:2],ddof=1)/2 + .75**2*np.var(values[2:],ddof=1)/3
    np.testing.assert_allclose(result['mc_standard_error']**2,variance)
