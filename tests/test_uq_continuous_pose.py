import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous import cell_forward, continuous_residual_norm
from fourier_splats.uq_continuous_pose import (PAIRS, PAIR_SCALE, cube_quadrature, polynomial_kernel_error,
    pose_derivative_fields, derivative_coefficient_bounds, grouped_gram_bounds, continuous_pose_audit,
    pose_cell_forward)


def test_continuous_pose_derivatives_against_nonlinear_fields():
    rng = np.random.default_rng(609390)
    k = rng.normal(size=(7, 3)); q = rng.normal(size=(7, 2)); c = rng.normal(size=7)+1j*rng.normal(size=7)
    xyz = rng.uniform(-.5, .5, (73, 3)); angle = .06; shift = .007
    d1, d2 = pose_derivative_fields(k, q, c, xyz, angle, shift, 'direct')
    fast1, fast2 = pose_derivative_fields(k, q, c, xyz, angle, shift, 'nufft')
    np.testing.assert_allclose(fast1, d1, atol=1e-10)
    np.testing.assert_allclose(fast2, d2, atol=1e-10)
    def field(u):
        rotated = k@Rotation.from_rotvec(angle*u[:3]).as_matrix()
        return (c@np.exp(2j*np.pi*(rotated@xyz.T+shift*(q@u[3:])[:, None]))).real
    u = rng.normal(size=5); u /= np.linalg.norm(u); step = .002
    symmetric = np.array([s*u[a]*u[b] for (a, b), s in zip(PAIRS, PAIR_SCALE)])
    np.testing.assert_allclose(d1@u, (field(step*u)-field(-step*u))/(2*step), atol=1e-7)
    np.testing.assert_allclose(d2@symmetric, (field(step*u)+field(-step*u)-2*field(np.zeros(5)))/step**2, atol=1e-7)
    assert np.isclose(np.linalg.norm(symmetric), 1.)


def test_polynomial_quadrature_bounds_and_grouped_fields():
    rng = np.random.default_rng(609391)
    k = rng.normal(size=(4, 3)); q = rng.normal(size=(4, 2)); c = rng.normal(size=4)+1j*rng.normal(size=4)
    angle = .1; shift = .01
    low_x, low_w = cube_quadrature(3); high_x, high_w = cube_quadrature(25)
    def design(x):
        d1, d2 = pose_derivative_fields(k, q, c, x, angle, shift, 'direct')
        return np.concatenate([d1, .5*d2], axis=1)
    low = design(low_x); high = design(high_x)
    low_gram = low.T@(low_w[:, None]*low); high_gram = high.T@(high_w[:, None]*high)
    first, second = derivative_coefficient_bounds(k, q, c, angle, shift)
    sums = np.r_[first, .5*second]; error = polynomial_kernel_error(k, 3)
    assert np.all(abs(low_gram-high_gram) <= error*sums[:, None]*sums[None, :]+1e-10)
    bounds = grouped_gram_bounds(low_gram, sums, [5, 15], error)
    for _ in range(30):
        u = rng.normal(size=5); u /= np.linalg.norm(u)
        v = np.r_[u, [s*u[a]*u[b] for (a, b), s in zip(PAIRS, PAIR_SCALE)]]
        assert np.sqrt(v@high_gram@v) <= bounds['minimum']+1e-10
    # High-order integration converges well before its explicit bound matters.
    higher_x, higher_w = cube_quadrature(30); higher = design(higher_x)
    np.testing.assert_allclose(higher.T@(higher_w[:, None]*higher), high_gram, atol=1e-11)


def test_continuous_pose_audit_covers_independent_cell_generators():
    rng = np.random.default_rng(609392); n, nq, box = 3, 4, 8
    k = rng.normal(size=(n, nq, 3)); q = rng.normal(size=(n, nq, 2)); ctf = rng.normal(size=(n, nq))
    noise = .8; weights = rng.normal(size=2*n*nq); angle = .15; shift = .02; B = .7
    pilot = rng.normal(size=box**3); pilot /= np.linalg.norm(pilot)
    residual = continuous_residual_norm(k, ctf, weights, noise, [[0, 0, 0]], [1], .17)
    audit = continuous_pose_audit(k, q, ctf, weights, noise, angle, shift, B, 1., residual['residual_norm'], order=24)
    assert audit['total_bias'] > audit['density_bias']
    nominal = cell_forward(k, ctf, pilot, box, noise)
    from fourier_splats.uq_continuous import cell_target_coefficients
    ell = cell_target_coefficients(box, [[0, 0, 0]], [1], .17)
    for _ in range(20):
        pose = rng.normal(size=(n, 5)); pose /= np.linalg.norm(pose, axis=1)[:, None]
        delta = rng.normal(size=box**3); delta *= B/np.linalg.norm(delta)
        observed = pose_cell_forward(k, q, ctf, pilot+delta, box, noise, pose, angle, shift)
        bias = weights@(observed-nominal)-ell@delta
        assert abs(bias) <= audit['total_bias']+1e-9
    zero = np.zeros((n, 5))
    np.testing.assert_allclose(pose_cell_forward(k, q, ctf, pilot, box, noise, zero, angle, shift), nominal, atol=1e-10)


def test_continuous_nonlinear_adversary_and_cubic_remainder():
    from fourier_splats.uq_continuous_pose import perturbed_geometry
    rng = np.random.default_rng(609394); n, nq, box = 2, 5, 8
    k = rng.normal(size=(n, nq, 3)); q = rng.normal(size=(n, nq, 2)); ctf = rng.normal(size=(n, nq))
    noise = .9; weights = rng.normal(size=2*n*nq); angle = .11; shift = .006; B = .8
    pilot = rng.normal(size=box**3); pilot /= np.linalg.norm(pilot)
    nominal = cell_forward(k, ctf, pilot, box, noise)
    residual = continuous_residual_norm(k, ctf, weights, noise, [[0, 0, 0]], [1], .2)
    audit = continuous_pose_audit(k, q, ctf, weights, noise, angle, shift, B, 1., residual['residual_norm'], order=24)
    xyz, quad = cube_quadrature(30)
    w = weights.reshape(n, 2*nq); coefficients = (w[:, :nq]+1j*w[:, nq:])*ctf/noise
    for _ in range(12):
        pose = rng.normal(size=(n, 5)); pose /= np.linalg.norm(pose, axis=1)[:, None]
        rotated, translations = perturbed_geometry(k, q, pose, angle, shift)
        moved = pose_cell_forward(k, q, ctf, pilot, box, noise, pose, angle, shift)
        pilot_bias = float(weights@(moved-nominal))
        exact = continuous_residual_norm(rotated, ctf*np.exp(2j*np.pi*translations), weights, noise, [[0, 0, 0]], [1], .2)
        # This eliminates all L2 density errors, including non-cell functions.
        assert abs(pilot_bias)+B*exact['residual_norm'] <= audit['total_bias']+1e-9
        error_field = np.zeros(len(xyz))
        for i in range(n):
            d1, d2 = pose_derivative_fields(k[i], q[i], coefficients[i], xyz, angle, shift, 'direct')
            sym = np.array([s*pose[i, a]*pose[i, b] for (a, b), s in zip(PAIRS, PAIR_SCALE)])
            actual = (coefficients[i]@np.exp(2j*np.pi*(rotated[i]@xyz.T+translations[i, :, None]))).real
            base = (coefficients[i]@np.exp(2j*np.pi*(k[i]@xyz.T))).real
            error_field += actual-base-d1@pose[i]-.5*d2@sym
        # Independent dense integration checks the uniform Taylor remainder.
        assert (1+B)*np.sqrt(quad@error_field**2) <= audit['remainder_bias']+1e-9
