import numpy as np
import pytest
from scipy.stats import norm
from fourier_splats.uq_continuous_gaussian import (
    continuous_gaussian, gaussian_ball_coverage, cell_pose_jacobian)
from fourier_splats.uq_fourier_pose_baseline import GaussianNuisanceWhitening
from fourier_splats.uq_continuous import ContinuousObservationGram
from fourier_splats.uq_continuous_pose import pose_cell_forward


class DenseGram:
    def __init__(self, a, ell):
        self.a = a
        self.ell = ell
        self.g = a@a.T
        self.shape = self.g.shape
        self.diagonal = np.diag(self.g)
    def target(self, *args):
        return self.a@self.ell, float(self.ell@self.ell)
    def matvec(self, w):
        return self.g@w


@pytest.mark.parametrize('tau', [.3, 1., 2.])
@pytest.mark.parametrize('pose', [False, True])
def test_matches_dense_conditional_gaussian(tau, pose):
    rng = np.random.default_rng(700)
    a = rng.normal(size=(24, 9))
    ell = rng.normal(size=9)
    nuisance = GaussianNuisanceWhitening(rng.normal(size=(3, 8, 5))*.3) if pose else None
    sigma = np.eye(24) if not pose else np.column_stack([nuisance.covariance_apply(x) for x in np.eye(24)])
    expected_w = np.linalg.solve(tau*tau*a@a.T+sigma, tau*tau*a@ell)
    # Independently condition in density space, using the precision identity.
    posterior = np.linalg.inv(np.eye(9)/tau**2+a.T@np.linalg.solve(sigma, a))
    fit = continuous_gaussian(DenseGram(a, ell), [[0, 0, 0]], [1], .1, tau, cg_rtol=1e-13)
    fit_pose = continuous_gaussian(DenseGram(a, ell), [[0, 0, 0]], [1], .1, tau, nuisance=nuisance, cg_rtol=1e-13)
    np.testing.assert_allclose(fit_pose['weights'], expected_w, rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(fit_pose['variational_variance_upper'], ell@posterior@ell, rtol=1e-9)
    assert fit_pose['converged']
    assert fit_pose['variance_identity_relative_error'] < 1e-8
    assert fit_pose['half_width'] >= fit['half_width']-1e-9


def test_nonconvergence_is_returned_not_hidden():
    rng = np.random.default_rng(2)
    fit = continuous_gaussian(DenseGram(rng.normal(size=(32, 21)), rng.normal(size=21)),
        [[0, 0, 0]], [1], .1, 2., maxiter=1)
    assert not fit['converged'] and fit['cg_info'] > 0
    assert np.isfinite(fit['half_width'])


def test_normal_envelope_uniform_coverage_and_threshold_failure():
    # Checks the analytic series proof over extreme ratios, including pure bias/noise.
    for alpha in [.05/12, .05, 2*norm.sf(np.sqrt(3))]:
        z = norm.isf(alpha/2)
        for t in np.r_[0., np.geomspace(1e-5, 1e5, 180)]:
            assert gaussian_ball_coverage(1., t, z*np.hypot(1., t)) >= 1-alpha-2e-14
    assert gaussian_ball_coverage(0., 2., 2.) == 1.
    assert gaussian_ball_coverage(1., .5, norm.isf(.25/2)*np.hypot(1., .5)) < .75


def test_continuous_pilot_jacobian_against_nonlinear_cells():
    rng = np.random.default_rng(712)
    k = rng.uniform(-2, 2, size=(2, 7, 3))
    q = rng.uniform(-2, 2, size=(2, 7, 2))
    ctf = rng.normal(size=(2, 7)); rho = rng.normal(size=8**3)
    j = cell_pose_jacobian(k, q, ctf, rho, 8, .7)
    h = 1e-5
    for a in range(5):
        u = np.zeros((2, 5)); u[:, a] = 1.
        plus = pose_cell_forward(k, q, ctf, rho, 8, .7, u, h, h)
        minus = pose_cell_forward(k, q, ctf, rho, 8, .7, -u, h, h)
        np.testing.assert_allclose(j[..., a], ((plus-minus)/(2*h)).reshape(2, 14), rtol=2e-7, atol=1e-8)


def test_exact_continuous_gram_and_ball_guarantee():
    rng = np.random.default_rng(91)
    gram = ContinuousObservationGram(rng.uniform(-2, 2, (2, 5, 3)), rng.normal(size=(2, 5)), .3)
    fit = continuous_gaussian(gram, [[.1, -.1, 0]], [1], .13, 2., density_radius=2.)
    assert fit['converged']
    assert fit['fixed_pose_uniform_class_coverage'] >= .95-1e-12
    for tau in [0., np.inf, np.nan, -1.]:
        with pytest.raises(ValueError):
            continuous_gaussian(gram, [[0, 0, 0]], [1], .1, tau, cg_rtol=1e-13)
