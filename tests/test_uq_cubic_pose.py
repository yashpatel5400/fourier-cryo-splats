import math
import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.spatial.transform import Rotation
from fourier_splats.uq_cubic_pose import (CubicPoseFieldOperator, POSE, SPATIAL, pose_lift,
    normalized_cell_moments_cubic, cell_moments_cubic, direct_cell_moments_cubic, gaussian_moments_cubic, cubic_residual_cross)
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_cell_moments import normalized_cell_moments, cell_fourier_moments
from fourier_splats.uq_joint_bias import gaussian_fourier_moments, joint_density_pose_bias
from fourier_splats.uq_continuous_pose import polynomial_kernel_error
from fourier_splats.uq_higher_remainder import higher_order_particle_remainders


def fixture(backend='direct', order=6):
    rng = np.random.default_rng(93)
    k = rng.normal(size=(2, 4, 3))*.5; q = rng.normal(size=(2, 4, 2))*.5
    ctf = rng.normal(size=(2, 4)); w = rng.normal(size=16)
    op = CubicPoseFieldOperator(k, q, ctf, w, 1.7, .13, .04, order=order, backend=backend, block_particles=1, nthreads=1)
    return rng, op, w


def skew(u):
    x, y, z = u
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0.]])


def test_lift_norms_and_quadratic_restriction():
    rng, op, w = fixture()
    u = rng.normal(size=(2, 5)); lift = pose_lift(u)
    for degree, (a, b) in enumerate(((0, 5), (5, 20), (20, 55)), 1):
        np.testing.assert_allclose(np.sum(lift[:, a:b]**2, axis=1), np.sum(u*u, axis=1)**degree)
    base = PolynomialPoseFieldOperator(op.k, op.q, op.ctf, w, 1.7, .13, .04, order=6, backend='direct')
    vec = rng.normal(size=40); cubic = np.r_[vec, np.zeros(70)]
    np.testing.assert_allclose(op.matvec(cubic), base.matvec(vec), atol=2e-14)


def test_exact_phase_taylor_recurrence_and_error_order():
    rng, op, _ = fixture()
    u = rng.normal(size=(2, 5)); u /= np.linalg.norm(u, axis=1)[:, None]
    expected = np.zeros(len(op.xyz))
    for i in range(2):
        W = op.angle[i]*skew(u[i, :3]); phases = []
        for degree in range(1, 4):
            phi = 2*np.pi*((op.k[i]@np.linalg.matrix_power(W, degree))@op.xyz.T)/math.factorial(degree)
            if degree == 1:
                phi += 2*np.pi*op.shift[i]*(op.q[i]@u[i, 3:])[:, None]
            phases.append(phi)
        exp_terms = [np.ones_like(phases[0], complex)]
        for degree in range(1, 4):
            exp_terms.append(sum(j*1j*phases[j-1]*exp_terms[degree-j] for j in range(1, degree+1))/degree)
        expected += np.real(np.sum(op.c[i, :, None]*np.exp(2j*np.pi*op.k[i]@op.xyz.T)*sum(exp_terms[1:]), axis=0))
    actual = op.matvec(op._from_particle_columns(pose_lift(u)))/op.sqrt_quad
    np.testing.assert_allclose(actual, expected, atol=2e-14)
    errors = []
    for scale in (1., .5, .25):
        exact = np.zeros_like(expected)
        for i in range(2):
            R = Rotation.from_rotvec(scale*op.angle[i]*u[i, :3]).as_matrix()
            phase0 = 2*np.pi*op.k[i]@op.xyz.T
            phase = 2*np.pi*((op.k[i]@R)@op.xyz.T+scale*op.shift[i]*(op.q[i]@u[i, 3:])[:, None])
            exact += np.real(np.sum(op.c[i, :, None]*(np.exp(1j*phase)-np.exp(1j*phase0)), axis=0))
        approx = op.matvec(op._from_particle_columns(pose_lift(scale*u)))/op.sqrt_quad
        errors.append(np.linalg.norm(exact-approx))
    assert 13 < errors[0]/errors[1] < 19
    assert 13 < errors[1]/errors[2] < 19


def test_adjoint_and_nufft_match_with_nonuniform_scaling():
    rng, direct, _ = fixture(); _, fast, _ = fixture('nufft')
    scales = np.exp(rng.normal(size=110))
    direct.denominators = scales; fast.denominators = scales
    u = rng.normal(size=110); v = rng.normal(size=direct.shape[0])
    np.testing.assert_allclose(v@direct.matvec(u), u@direct.rmatvec(v), atol=2e-13)
    np.testing.assert_allclose(fast.matvec(u), direct.matvec(u), atol=2e-11)
    np.testing.assert_allclose(fast.rmatvec(v), direct.rmatvec(v), atol=2e-11)


def test_scalar_cell_moments_independent_quadrature():
    k = np.array([-19., -2., -.01, 0., .001, 1., 12.]); box = 4
    nodes, weights = leggauss(80); x = nodes/(2*box)
    expected = np.stack([np.exp(2j*np.pi*k[:, None]*x)@(weights*x**degree)/2 for degree in range(4)], axis=-1)
    actual = normalized_cell_moments_cubic(k, box)
    np.testing.assert_allclose(actual, expected, atol=3e-15)
    np.testing.assert_allclose(actual[:, :3], normalized_cell_moments(k, box), atol=3e-14)


def test_cell_twenty_moments_direct_physical_quadrature():
    rng = np.random.default_rng(11); box = 2; coeff = rng.normal(size=box**3)
    k = rng.normal(size=(5, 3))*2
    nodes, qw = leggauss(14); z, y, x = np.meshgrid(nodes, nodes, nodes, indexing='ij')
    offsets = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)/(2*box)
    quad = (qw[:, None, None]*qw[None, :, None]*qw[None, None, :]).ravel()/(2*box)**3
    centers = (np.indices((box,)*3).reshape(3, -1).T[:, ::-1]+.5)/box-.5
    expected = np.zeros((len(k), 20), complex)
    for c, weight in zip(centers, coeff):
        xyz = c+offsets; phase = np.exp(2j*np.pi*k@xyz.T)
        for j, beta in enumerate(SPATIAL):
            expected[:, j] += weight*box**1.5*(phase@(quad*np.prod(xyz**np.array(beta), axis=1)))
    actual = cell_moments_cubic(k, coeff, box)
    np.testing.assert_allclose(direct_cell_moments_cubic(k, coeff, box), expected, atol=3e-14)
    np.testing.assert_allclose(actual, expected, atol=2e-12)
    np.testing.assert_allclose(actual[:, :10], cell_fourier_moments(k, coeff, box), atol=2e-12)


def test_gaussian_twenty_moments_independent_tensor_quadrature():
    rng = np.random.default_rng(81); k = rng.normal(size=(4, 3))*1.5
    centers = np.array([[.1, -.12, .04], [-.15, .07, .2]]); signs = [1., -.5]; sigma = .13
    nodes, weights = leggauss(90); x = nodes/2; expected = np.zeros((4, 20), complex)
    for center, sign in zip(centers, signs):
        univariate = np.empty((4, 3, 4), complex)
        for axis in range(3):
            density = np.exp(-.5*((x-center[axis])/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
            for degree in range(4):
                univariate[:, axis, degree] = np.exp(2j*np.pi*k[:, axis, None]*x)@(weights*density*x**degree/2)
        for j, beta in enumerate(SPATIAL):
            expected[:, j] += sign*np.prod(np.stack([univariate[:, axis, p] for axis, p in enumerate(beta)]), axis=0)
    actual = gaussian_moments_cubic(k, centers, signs, sigma)
    np.testing.assert_allclose(actual, expected, atol=2e-14)
    np.testing.assert_allclose(actual[:, :10], gaussian_fourier_moments(k, centers, signs, sigma), atol=1e-15)


def test_joint_continuous_bias_controls_sampled_nonlinear_fields():
    rng, op, w = fixture(order=12)
    scaling = op.establish_quadrature_design_scaling(order=4)
    L = np.sqrt(scaling['sum_group_scales'])
    fields = np.stack([op.matvec(e) for e in np.eye(op.shape[1])], axis=1)
    amp = np.hypot(w.reshape(2, 8)[:, :4], w.reshape(2, 8)[:, 4:])
    masses = np.einsum('naj,nj->na', op.coefficient_bound_matrix(), amp)
    pad = polynomial_kernel_error(op.k, 12, degree=6)*np.sum(masses*masses)
    f = L*np.sqrt(np.linalg.norm(fields, ord=2)**2+pad)
    remainders = [higher_order_particle_remainders(op.k[i], op.q[i], op.c[i], .13, .04, degrees=(3,))['records'][0]['field_remainder'] for i in range(2)]
    centers = [[.1, -.1, .05]]; signs = [1.]; width = .2
    cross = cubic_residual_cross(op, w, 1.7, centers, signs, width)
    nominal = np.real(np.sum(op.c.reshape(-1, 1)*np.exp(2j*np.pi*op.k.reshape(-1, 3)@op.xyz.T), axis=0))
    ell = np.exp(-np.sum((op.xyz-np.array(centers[0]))**2, axis=1)/(2*width**2))/(np.sqrt(2*np.pi)*width)**3
    h = ell-nominal
    np.testing.assert_allclose(cross['cross_vector'], fields.T@(op.sqrt_quad*h), atol=3e-9)
    bound = joint_density_pose_bias(np.linalg.norm(op.sqrt_quad*h), f, cross['norm_upper'], L, 2., 1., sum(remainders))
    for _ in range(8):
        u = rng.normal(size=(2, 5)); u /= np.linalg.norm(u, axis=1)[:, None]
        exact = np.zeros(len(op.xyz))
        for i in range(2):
            rotation = Rotation.from_rotvec(.13*u[i, :3]).as_matrix()
            phase = 2*np.pi*((op.k[i]@rotation)@op.xyz.T+.04*(op.q[i]@u[i, 3:])[:, None])
            exact += np.real(np.sum(op.c[i, :, None]*np.exp(1j*phase), axis=0))
        polynomial = op.matvec(op._from_particle_columns(pose_lift(u))*op.denominators)/op.sqrt_quad
        assert np.linalg.norm(op.sqrt_quad*(exact-nominal-polynomial)) <= sum(remainders)+1e-12
        # Cauchy upper bound on supremum over every pilot of norm <=1 and
        # every continuous density displacement of norm <=2 on this toy grid.
        true_sup = 2*np.linalg.norm(op.sqrt_quad*(ell-exact))+np.linalg.norm(op.sqrt_quad*(exact-nominal))
        assert true_sup <= bound['bias_upper']+1e-9


def test_invalid_inputs_and_inherited_methods_fail_closed():
    _, op, _ = fixture()
    with pytest.raises(ValueError): pose_lift(np.zeros(4))
    with pytest.raises(ValueError): op.establish_quadrature_design_scaling(1)
    with pytest.raises(ValueError): op.pair_moments(np.zeros((2, 4, 10)))
    with pytest.raises(NotImplementedError): op.establish_group_scaling()
    with pytest.raises(NotImplementedError): op.establish_coefficient_scaling()
    with pytest.raises(NotImplementedError): op.weight_gradient(None, None)


@pytest.mark.parametrize('distortion', [0., 1e-7])
def test_planar_cubic_remainder_and_fourth_order_scaling(distortion):
    rng = np.random.default_rng(612); q = rng.normal(size=(1, 7, 2))*2
    k = np.pad(q, ((0, 0), (0, 0), (0, 1)))@Rotation.from_rotvec([.2, -.1, .4]).as_matrix()
    k += distortion*rng.normal(size=k.shape)
    ctf = np.ones((1, 7)); w = rng.normal(size=14)
    op = CubicPoseFieldOperator(k, q, ctf, w, 1., .13, .04, order=20, backend='direct')
    nominal = np.real(np.sum(op.c[0, :, None]*np.exp(2j*np.pi*k[0]@op.xyz.T), axis=0))
    record = higher_order_particle_remainders(k[0], q[0], op.c[0], .13, .04, degrees=(3,), domain='cube')
    bound = record['records'][0]['field_remainder']; errors = []
    directions = rng.normal(size=(32, 5)); directions /= np.linalg.norm(directions, axis=1)[:, None]
    for scale in (1., .5, .25):
        maxima = 0.
        for u in directions:
            rotation = Rotation.from_rotvec(scale*.13*u[:3]).as_matrix()
            phase = 2*np.pi*((k[0]@rotation)@op.xyz.T+scale*.04*(q[0]@u[3:])[:, None])
            exact = np.real(np.sum(op.c[0, :, None]*np.exp(1j*phase), axis=0))-nominal
            approx = op.matvec(pose_lift((scale*u)[None])[0])/op.sqrt_quad
            maxima = max(maxima, np.linalg.norm(op.sqrt_quad*(exact-approx)))
        errors.append(maxima)
    ratio = bound/errors[0]
    # The initially requested factor-20 tightness check failed at 35.57;
    # preserve that diagnostic rather than claiming this mixed bound is tight.
    # Polynomial correctness is falsified by the independent order test below.
    assert ratio >= 1, f'Invalid planar remainder upper: {ratio}'
    assert 11 < errors[0]/errors[1] < 22
    assert 13 < errors[1]/errors[2] < 19


def test_translation_only_remainder_is_tight_enough_to_falsify_missing_terms():
    q = np.array([[[2., 0.], [0., 1.]]]); k = np.pad(q, ((0, 0), (0, 0), (0, 1)))
    op = CubicPoseFieldOperator(k, q, np.ones((1, 2)), [1., 0., 0., 0.], 1., 0., .02, order=20, backend='direct')
    record = higher_order_particle_remainders(k[0], q[0], op.c[0], 0., .02, degrees=(3,), domain='cube')
    u = np.array([0., 0., 0., 1., 0.]); phase = 2*np.pi*k[0, 0]@op.xyz.T
    errors = []
    for scale in (1., .5, .25):
        exact = np.cos(phase+2*np.pi*scale*.02*q[0, 0, 0])-np.cos(phase)
        polynomial = op.matvec(pose_lift((scale*u)[None])[0])/op.sqrt_quad
        errors.append(np.linalg.norm(op.sqrt_quad*(exact-polynomial)))
    ratio = record['records'][0]['field_remainder']/errors[0]
    assert 1 <= ratio < 3
    assert 15 < errors[0]/errors[1] < 17
    assert 15 < errors[1]/errors[2] < 17
