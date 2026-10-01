import numpy as np
from fourier_splats.uq_mixed_pose import first_derivative_gram, mixed_pose_interval
from fourier_splats.uq_continuous_pose import cube_quadrature, pose_derivative_fields, pose_cell_forward
from fourier_splats.uq_continuous import cell_forward, continuous_residual_norm
from fourier_splats.uq_continuous_gaussian import cell_pose_jacobian


def test_analytic_gram_against_independent_dense_quadrature():
    rng = np.random.default_rng(511)
    k = rng.normal(size=(8, 3)); q = rng.normal(size=(8, 2))
    c = rng.normal(size=8)+1j*rng.normal(size=8)
    xyz, qw = cube_quadrature(25)
    d, _ = pose_derivative_fields(k, q, c, xyz, .08, .02, 'direct')
    exact = d.T@(qw[:, None]*d)
    fit = first_derivative_gram(k, q, c, .08, .02)
    np.testing.assert_allclose(fit['gram'], exact, rtol=2e-11, atol=1e-11)
    assert fit['spectral_norm_upper'] >= np.sqrt(np.linalg.eigvalsh(exact)[-1])-1e-11


def test_independent_scaling_and_common_mode_floor():
    rng = np.random.default_rng(512)
    k = rng.normal(size=(1, 6, 3)); q = rng.normal(size=(1, 6, 2)); ctf = np.ones((1, 6))
    w = rng.normal(size=12)
    one = mixed_pose_interval(k, q, ctf, w, .3, .03, .01, 2., 1., .1, common_radius=.1)
    many = mixed_pose_interval(np.repeat(k, 9, axis=0), np.repeat(q, 9, axis=0),
        np.repeat(ctf, 9, axis=0), np.tile(w/9, 9), .3, .03, .01, 2., 1., .1, common_radius=.1)
    np.testing.assert_allclose(many['stochastic_pose_variance_proxy'], one['stochastic_pose_variance_proxy']/9)
    np.testing.assert_allclose(many['measurement_variance'], one['measurement_variance']/9)
    np.testing.assert_allclose(many['common_mode_bias'], one['common_mode_bias'])
    np.testing.assert_allclose(many['nonlinear_remainder_bias'], one['nonlinear_remainder_bias'])


def test_first_order_and_remainder_cover_nonlinear_cell_examples():
    rng = np.random.default_rng(513); n, nq, box = 3, 7, 8
    k = rng.normal(size=(n, nq, 3)); q = rng.normal(size=(n, nq, 2)); ctf = rng.normal(size=(n, nq))
    w = rng.normal(size=2*n*nq); noise = .7; angle = .08; shift = .02
    pilot = rng.normal(size=box**3); pilot /= np.linalg.norm(pilot)
    delta = rng.normal(size=box**3); delta *= .5/np.linalg.norm(delta); rho = pilot+delta
    residual = continuous_residual_norm(k, ctf, w, noise, [[0, 0, 0]], [1], .2)['residual_norm']
    j = cell_pose_jacobian(k, q, ctf, pilot, box, noise)*np.array([angle]*3+[shift]*2)
    pair = np.einsum('nmp,nm->np', j, w.reshape(n, -1))
    result = mixed_pose_interval(k, q, ctf, w, noise, angle, shift, .5, 1., residual,
        independent_radius=.7, common_radius=.2, pilot_derivative_pairings=pair)
    jr = cell_pose_jacobian(k, q, ctf, rho, box, noise)*np.array([angle]*3+[shift]*2)
    true_pair = np.einsum('nmp,nm->np', jr, w.reshape(n, -1))
    assert np.all(np.linalg.norm(true_pair, axis=1) <= np.array(result['scalar_derivative_bounds'])+1e-10)
    nominal = cell_forward(k, ctf, rho, box, noise)
    for _ in range(10):
        u = rng.normal(size=(n, 5)); u *= .9/np.linalg.norm(u, axis=1)[:, None]
        nonlinear = float(w@(pose_cell_forward(k, q, ctf, rho, box, noise, u, angle, shift)-nominal))
        linear = float(np.sum(true_pair*u))
        assert abs(nonlinear-linear) <= result['nonlinear_remainder_bias']+1e-10
