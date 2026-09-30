import numpy as np
from numpy.polynomial.legendre import leggauss
from fourier_splats.uq_polynomial_trace import interval_polynomial_moments, polynomial_fourier_block_trace
from fourier_splats.uq_cubic_pose import SPATIAL
from fourier_splats.uq_continuous_pose import cube_quadrature


def test_seventh_scalar_moments_against_oscillatory_independent_quadrature():
    frequencies = np.array([-24., -12., -.2, 0., 1e-10, 3., 24.]); x, w = leggauss(240); x /= 2; w /= 2
    actual = interval_polynomial_moments(frequencies)
    expected = np.stack([np.exp(2j*np.pi*frequencies[:, None]*x)@(w*x**j) for j in range(7)], axis=-1)
    # At integer frequencies the exact zeroth moment is zero; the 240-node
    # quadrature has 5.8e-15 absolute cancellation error at +/-24.
    np.testing.assert_allclose(actual, expected, atol=8e-15)
    assert np.max(abs(actual[[0, 6], 0])) < 4e-16


def test_trace_matches_independent_spatial_quadrature_with_mixed_complex_coefficients():
    rng = np.random.default_rng(611); k = rng.normal(size=(5, 3))*1.3
    c = rng.normal(size=(5, 3, 20))+1j*rng.normal(size=(5, 3, 20))
    xyz, q = cube_quadrature(24); monomials = np.stack([np.prod(xyz**np.array(beta), axis=1) for beta in SPATIAL])
    phase = np.exp(2j*np.pi*k@xyz.T)
    values = np.einsum('kap,kr,pr->ra', c, phase, monomials).real
    expected = np.sum(q[:, None]*values**2)
    result = polynomial_fourier_block_trace(k, c, SPATIAL)
    np.testing.assert_allclose(result['trace_unpadded'], expected, rtol=3e-13, atol=1e-12)
    assert result['heuristic_roundoff_guard'] >= 0
