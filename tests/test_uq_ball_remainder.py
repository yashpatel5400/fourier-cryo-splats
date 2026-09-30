import itertools
import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.spatial.transform import Rotation
from fourier_splats.uq_ball_remainder import ball_fourier_kernel, cube_fourier_kernel, ball_sobolev_norms, particle_ball_remainder
from fourier_splats.uq_continuous_pose import cube_quadrature, pose_derivative_fields, PAIRS, PAIR_SCALE


def spherical_quadrature(radius):
    r, rw = leggauss(28); r = (r+1)*radius/2; rw = rw*radius/2*r*r
    z, zw = leggauss(36); phi = np.arange(80)*2*np.pi/80
    rr, zz, pp = np.meshgrid(r, z, phi, indexing='ij')
    xyz = np.stack([rr*np.sqrt(1-zz*zz)*np.cos(pp), rr*np.sqrt(1-zz*zz)*np.sin(pp), rr*zz], -1).reshape(-1, 3)
    weight = np.broadcast_to(rw[:, None, None]*zw[None, :, None]*2*np.pi/80, rr.shape).ravel()
    return xyz, weight


def test_ball_kernel_and_all_ordered_derivatives_against_spherical_integration():
    rng = np.random.default_rng(610401); k = rng.normal(size=(7, 3))*.55
    c = rng.normal(size=7)+1j*rng.normal(size=7); radius = .91
    xyz, weight = spherical_quadrature(radius)
    phase = np.exp(2j*np.pi*(xyz@k.T))
    np.testing.assert_allclose(ball_fourier_kernel(k, radius), weight@phase, atol=3e-14)
    assert ball_fourier_kernel(np.zeros((1, 3)), radius)[0] == pytest.approx(4*np.pi*radius**3/3)
    norms, diagnostics = ball_sobolev_norms(k, c, radius, orders=(0, 1, 2, 3))
    for m, norm, diagnostic in zip(range(4), norms, diagnostics):
        squared = 0.
        for axes in itertools.product(range(3), repeat=m):
            coefficient = c*np.prod((2j*np.pi*k[:, axes]), axis=1) if axes else c
            field = (phase@coefficient).real
            squared += weight@(field*field)
        np.testing.assert_allclose(norm**2-diagnostic['roundoff_pad'], squared, rtol=2e-12, atol=2e-12)


@pytest.mark.parametrize('distortion', [0., .03])
@pytest.mark.parametrize('domain', ['ball', 'cube'])
def test_nonlinear_cube_remainder_with_joint_poses_and_embedding_residual(distortion, domain):
    rng = np.random.default_rng(610403); q = rng.normal(size=(9, 2))
    embedding = Rotation.from_rotvec([.2, -.3, .1]).as_matrix()[:2]
    k = q@embedding+distortion*rng.normal(size=(9, 3))
    c = rng.normal(size=9)+1j*rng.normal(size=9)
    angle, shift = .17, .035
    bound = particle_ball_remainder(k, q, c, angle, shift, domain=domain)
    xyz, quadrature = cube_quadrature(26)
    nominal = (np.exp(2j*np.pi*(xyz@k.T))@c).real
    first, second = pose_derivative_fields(k, q, c, xyz, angle, shift, backend='direct')
    for _ in range(8):
        u = rng.normal(size=5); u /= np.linalg.norm(u)
        rotated = k@Rotation.from_rotvec(angle*u[:3]).as_matrix()
        actual = (np.exp(2j*np.pi*(xyz@rotated.T+shift*(q@u[3:])[None, :]))@c).real
        products = np.array([u[a]*u[b]*s for (a, b), s in zip(PAIRS, PAIR_SCALE)])
        remainder = actual-nominal-first@u-.5*second@products
        actual_norm = np.sqrt(quadrature@(remainder*remainder))
        assert actual_norm <= bound['field_remainder']*(1+1e-12)
    if distortion:
        assert bound['embedding_residual_remainder'] > 0


def test_expanded_cube_derivative_norms_against_direct_quadrature():
    rng = np.random.default_rng(610405); k = rng.normal(size=(8, 3))*.9
    c = rng.normal(size=8)+1j*rng.normal(size=8); radius = .527
    xyz, weights = cube_quadrature(32); xyz *= 2*radius; weights *= (2*radius)**3
    phase = np.exp(2j*np.pi*(xyz@k.T))
    np.testing.assert_allclose(cube_fourier_kernel(k, radius), weights@phase, atol=5e-14)
    norms, diagnostics = ball_sobolev_norms(k, c, radius, domain='cube')
    for m, norm, diagnostic in zip(range(1, 4), norms, diagnostics):
        squared = 0.
        for axes in itertools.product(range(3), repeat=m):
            field = (phase@(c*np.prod(2j*np.pi*k[:, axes], axis=1))).real
            squared += weights@(field*field)
        np.testing.assert_allclose(norm**2-diagnostic['roundoff_pad'], squared, rtol=3e-12)


def test_zero_pose_and_input_validation():
    q = np.array([[1., 0], [0, 1], [1, 1]])
    k = np.column_stack([q, np.zeros(3)])
    assert particle_ball_remainder(k, q, np.ones(3), 0., 0.)['field_remainder'] == 0.
    with pytest.raises(ValueError):
        ball_sobolev_norms(k, np.array([1., np.nan, 1.]), 1.)
    with pytest.raises(ValueError):
        particle_ball_remainder(k, q*0, np.ones(3), .1, .1)


def test_near_tight_pure_translation_against_closed_form_remainder():
    q = np.array([[.1, 0], [0, .2], [.1, .2]])
    k = np.column_stack([q, np.zeros(3)]); c = np.array([-1j, 0, 0])
    shift = .01; phase = 2*np.pi*q[0, 0]*shift
    residual_coefficient = c[0]*(np.expm1(1j*phase)-1j*phase+phase**2/2)
    actual_squared = .5*(abs(residual_coefficient)**2+(residual_coefficient**2*np.prod(np.sinc(2*k[0]))).real)
    bound = particle_ball_remainder(k, q, c, 0., shift, domain='cube')['field_remainder']
    assert 1. <= bound/np.sqrt(actual_squared) < 1.04


def test_antipodal_frequency_sign_against_closed_form_cosine_norms():
    k = np.array([[.07, -.1, .03], [-.07, .1, -.03]]); c = np.array([.5, .5]); radius = .53
    norms, diagnostics = ball_sobolev_norms(k, c, radius, orders=(np.int64(1), 2, 3), domain='cube')
    for m, value, diag in zip(range(1, 4), norms, diagnostics):
        exact = .5*(2*np.pi*np.linalg.norm(k[0]))**(2*m)*((2*radius)**3+(-1)**m*cube_fourier_kernel(2*k[:1], radius)[0])
        np.testing.assert_allclose(value**2-diag['roundoff_pad'], exact, rtol=1e-13)


@pytest.mark.parametrize('angle', [.2, np.pi, 1.2*np.pi])
def test_scaled_embedding_product_set_and_deterministic_pose_directions(angle):
    q = np.array([[.4, 0], [0, .5], [.3, -.2]])
    k = np.column_stack([.3*q, np.zeros(3)]); c = np.array([1.+.7j, -.8+.1j, .3-.4j]); shift = .08
    joint = particle_ball_remainder(k, q, c, angle, shift, domain='cube')
    product = particle_ball_remainder(k, q, c, angle, shift, joint_ball=False, domain='cube')
    np.testing.assert_allclose(product['translation_radius'], shift/.3)
    assert joint['path_speed_bound'] < product['path_speed_bound']
    assert joint['field_remainder'] < product['field_remainder']
    xyz, weights = cube_quadrature(18); phase = np.exp(2j*np.pi*(xyz@k.T)); nominal = (phase@c).real
    for axis in np.eye(3):
        for theta in np.arange(8)*np.pi/4:
            shift_direction = np.array([np.cos(theta), np.sin(theta)])
            cross = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])*angle
            p1 = 2*np.pi*(xyz@(k@cross).T+shift*(q@shift_direction)[None, :])
            p2 = 2*np.pi*(xyz@(k@cross@cross).T)
            first = (phase*(1j*p1)@c).real
            second = (phase*(1j*p2-p1*p1)@c).real
            rotated = k@Rotation.from_rotvec(angle*axis).as_matrix()
            actual = (np.exp(2j*np.pi*(xyz@rotated.T+shift*(q@shift_direction)[None, :]))@c).real
            error = actual-nominal-first-.5*second
            assert np.sqrt(weights@(error*error)) <= product['field_remainder']
            transformed = xyz@Rotation.from_rotvec(angle*axis).as_matrix().T+np.r_[shift*shift_direction/.3, 0]
            assert np.max(abs(transformed)) <= product['domain_radius']


def test_derivative_fields_against_independent_finite_differences():
    rng = np.random.default_rng(610407); q = rng.normal(size=(7, 2)); k = rng.normal(size=(7, 3))*.4
    c = rng.normal(size=7)+1j*rng.normal(size=7); xyz = rng.uniform(-.5, .5, size=(19, 3))
    angle, shift = .3, .1; u = rng.normal(size=5); u /= np.linalg.norm(u)
    def field(t):
        rotated = k@Rotation.from_rotvec(t*angle*u[:3]).as_matrix()
        return (np.exp(2j*np.pi*(xyz@rotated.T+t*shift*(q@u[3:])[None, :]))@c).real
    first, second = pose_derivative_fields(k, q, c, xyz, angle, shift, backend='direct')
    products = np.array([u[a]*u[b]*s for (a, b), s in zip(PAIRS, PAIR_SCALE)])
    h = .0002
    np.testing.assert_allclose(first@u, (field(h)-field(-h))/(2*h), rtol=1e-6, atol=2e-8)
    np.testing.assert_allclose(second@products, (field(h)-2*field(0)+field(-h))/h**2, rtol=1e-5, atol=2e-7)


def test_source_class_guard_rejects_product_or_unrecorded_density_class():
    import importlib.util
    from pathlib import Path
    spec = importlib.util.spec_from_file_location('enclosure_probe', Path(__file__).resolve().parents[1]/'scripts/probe_uq_ball_remainder.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    joint = {'pose_set': 'Per-particle joint ball: squared scaled rotation norm plus squared scaled translation norm <= 1',
             'density_radius': 2., 'pilot_norm_bound': 1.}
    assert module.validate_source_class(joint, {})[:2] == (2., 1.)
    with pytest.raises(ValueError):
        module.validate_source_class(dict(joint, pose_set='product'), {})
    with pytest.raises(ValueError):
        module.validate_source_class({'pose_set': joint['pose_set']}, {})
