"""Independent integration checks of uniform continuous remainder bounds."""
import numpy as np
from scipy.linalg import expm
from fourier_splats.uq_continuous_pose import cube_quadrature
from fourier_splats.uq_continuous_moments import continuous_third_derivative_bound, integrated_cubic_remainder


def test_cube_linear_moments_against_gaussian_envelope():
    xyz, weights = cube_quadrature(5)
    vectors = np.vstack([np.eye(3), np.ones((1, 3))/np.sqrt(3),
                         np.random.default_rng(321).normal(size=(40, 3))])
    for p, moment in [(2, 1), (4, 3), (6, 15)]:
        integral = np.sum(weights[:, None]*(xyz@vectors.T)**p, axis=0)
        envelope = moment*(np.linalg.norm(vectors, axis=1)**2/12)**(p/2)
        assert np.all(integral <= envelope+1e-12)
    np.testing.assert_allclose(np.sum(weights*xyz[:, 0]**6), 1/448, rtol=1e-12)


def test_integrated_bound_covers_exact_rotation_third_derivatives():
    xyz, weights = cube_quadrature(8)
    rng = np.random.default_rng(322)
    k = rng.normal(size=(24, 3))*5; q = rng.normal(size=(24, 2))*5
    for angle, shift in [(0., .02), (.02, 0.), (.12, .03), (.3, .05)]:
        envelope = continuous_third_derivative_bound(k, q, angle, shift)
        for _ in range(12):
            u = rng.normal(size=5); u /= np.linalg.norm(u)
            a, b, c = angle*u[:3]
            skew = np.array([[0., -c, b], [c, 0., -a], [-b, a, 0.]])
            rotation = expm(rng.uniform()*skew)
            first = 2*np.pi*((k@rotation@skew)@xyz.T+shift*(q@u[3:])[:, None])
            second = 2*np.pi*((k@rotation@skew@skew)@xyz.T)
            third = 2*np.pi*((k@rotation@skew@skew@skew)@xyz.T)
            # Unit-modulus exp(i phi) drops out of the squared derivative norm.
            derivative = 1j*third-3*first*second-1j*first**3
            integrated = np.sqrt(np.sum(weights[None]*abs(derivative)**2, axis=1))
            assert np.all(integrated <= envelope+1e-10)


def test_integrated_cubic_bound_dominates_exact_nonlinear_remainder():
    xyz, weights = cube_quadrature(32)
    rng = np.random.default_rng(323)
    n, nq = 3, 4
    k = rng.normal(size=(n, nq, 3))*3; q = rng.normal(size=(n, nq, 2))*3
    transfer = rng.uniform(-1, 1, (n, nq)); w = rng.normal(size=(n, 2*nq))
    noise = .7; angle = .12; shift = .025
    cert = integrated_cubic_remainder(k, q, transfer, w, noise, angle, shift, 2., 1.)
    residual = np.zeros(len(xyz))
    for i in range(n):
        u = rng.normal(size=5); u /= np.linalg.norm(u)
        a, b, c = angle*u[:3]
        skew = np.array([[0., -c, b], [c, 0., -a], [-b, a, 0.]])
        nominal_phase = 2*np.pi*(k[i]@xyz.T)
        phi1 = 2*np.pi*((k[i]@skew)@xyz.T+shift*(q[i]@u[3:])[:, None])
        phi2 = 2*np.pi*((k[i]@skew@skew)@xyz.T)
        final_phase = 2*np.pi*((k[i]@expm(skew))@xyz.T+shift*(q[i]@u[3:])[:, None])
        remainder = np.exp(1j*final_phase)-np.exp(1j*nominal_phase)*(1+1j*phi1+.5*(1j*phi2-phi1**2))
        coeff = (w[i, :nq]+1j*w[i, nq:])*transfer[i]/noise
        residual += (coeff@remainder).real
    actual_field_norm = np.sqrt(np.sum(weights*residual**2))
    assert 3*actual_field_norm <= cert['remainder_bias']+1e-10
    assert cert['frequency_weighted_remainder_bias'] <= cert['hilbert_schmidt_remainder_bias']+1e-10
    # Compare against the older pointwise spatial maximum as an ablation.
    radius = np.sqrt(3)/2
    L = 2*np.pi*(angle*np.linalg.norm(k, axis=-1)*radius+shift*np.linalg.norm(q, axis=-1))
    H = 2*np.pi*angle**2*np.linalg.norm(k, axis=-1)*radius
    T = 2*np.pi*angle**3*np.linalg.norm(k, axis=-1)*radius
    old = 3/6*np.dot(np.sqrt(np.sum((transfer/noise)**2*(L**3+3*L*H+T)**2, axis=1)),
                    np.linalg.norm(w, axis=1))
    assert cert['remainder_bias'] < old
