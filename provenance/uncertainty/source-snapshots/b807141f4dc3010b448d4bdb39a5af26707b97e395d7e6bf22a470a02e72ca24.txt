import numpy as np
import pytest
from scipy.special import logsumexp
from scipy.stats import norm
from fourier_splats.uq_mixture_validation import (
    gaussian_ball_log_upper, mixture_envelope_upper, log_validation_evalue)


def test_gaussian_ball_and_scale_envelope_against_independent_search():
    y = np.array([1.7, -.8, .2]); center = np.array([.3, .4, -.1])
    radius = .6; distance = np.linalg.norm(y-center)
    upper = gaussian_ball_log_upper(distance, 3, radius, .2, 2.)
    rng = np.random.default_rng(842)
    means = center + radius*rng.normal(size=(1000, 3)) / np.sqrt(3)
    means = center + (means-center)/np.maximum(np.linalg.norm(means-center, axis=1)[:, None]/radius, 1)
    scales = np.geomspace(.2, 2., 71)
    direct = norm.logpdf(y[None, None, :], loc=means[:, None, :], scale=scales[None, :, None]).sum(axis=-1)
    assert np.max(direct) <= upper+1e-12
    best_mean = center+radius*(y-center)/distance
    best_sigma = np.clip((distance-radius)/np.sqrt(3), .2, 2.)
    assert upper == pytest.approx(norm.logpdf(y, loc=best_mean, scale=best_sigma).sum())
    assert np.isposinf(gaussian_ball_log_upper(.2, 3, .3, 0., np.inf))
    assert log_validation_evalue(-3., np.inf) == -np.inf


def test_mixture_upper_against_independent_conic_optimizer_and_anchor():
    import cvxpy as cp
    rng = np.random.default_rng(643)
    kernels = np.exp(rng.normal(size=(19, 8)))
    fit = mixture_envelope_upper(np.log(kernels), tolerance=2e-6, max_iterations=10000)
    w = cp.Variable(8)
    problem = cp.Problem(cp.Maximize(cp.sum(cp.log(kernels@w))), [w >= 0, cp.sum(w) == 1])
    problem.solve(solver='CLARABEL', tol_gap_abs=1e-9, tol_gap_rel=1e-9, tol_feas=1e-9)
    assert problem.status == cp.OPTIMAL
    assert fit.converged
    assert fit.envelope_primal <= problem.value+1e-6
    assert fit.log_likelihood_upper >= problem.value-1e-6
    assert fit.log_likelihood_upper-problem.value < 5e-6
    anchor = kernels@fit.upper_anchor_weights
    replay = np.log(anchor).sum()-19+np.max((kernels/anchor[:, None]).sum(axis=0))
    assert fit.log_likelihood_upper == pytest.approx(replay, abs=1e-11)
    assert fit.weights.sum() == pytest.approx(1)


def test_early_stop_still_upper_bounds_all_discrete_mixtures():
    log_kernel = np.array([[0., -2., -4.], [-3., -.1, -7.], [-1., -2., -.1]])
    fit = mixture_envelope_upper(log_kernel, max_iterations=0)
    assert not fit.converged and fit.status == 'iteration_limit'
    for a in np.linspace(0, 1, 41):
        for b in np.linspace(0, 1-a, 31):
            w = np.array([a, b, 1-a-b])
            likelihood = np.log(np.exp(log_kernel)@w).sum()
            assert likelihood <= fit.log_likelihood_upper+1e-12
    shifted = mixture_envelope_upper(log_kernel+np.array([1000., -2000., -3000.])[:, None], max_iterations=0)
    assert shifted.log_likelihood_upper == pytest.approx(fit.log_likelihood_upper-4000.)
    assert np.isposinf(mixture_envelope_upper(np.array([[0., np.inf]])).log_likelihood_upper)
