import math
import numpy as np
import pytest
from scipy.linalg import expm
from fourier_splats.uq_continuous_pose import cube_quadrature
from fourier_splats.uq_higher_remainder import partial_bell, higher_order_particle_remainders


def test_fourth_order_chain_rule_coefficients_and_phase_bound():
    x = np.array([2., 3., 5., 7.]); b = partial_bell(x)
    np.testing.assert_allclose(b[4, 1:5], [7., 4*2*5+3*3**2, 6*2**2*3, 2**4])
    assert b[0, 0] == 1
    np.testing.assert_allclose(partial_bell(np.stack([x, 2*x], axis=1))[:, :, 0], b)


@pytest.mark.parametrize('nonplanar', [False, True])
@pytest.mark.parametrize('domain', ['ball', 'cube'])
def test_degree_two_matches_and_all_degrees_bound_direct_nonlinear_remainder(nonplanar, domain):
    rng = np.random.default_rng(609930)
    q = rng.normal(size=(7, 2))*.7
    k = np.column_stack([q, np.zeros(len(q))])
    if nonplanar:
        k += rng.normal(size=k.shape)*.035
    c = rng.normal(size=len(q))+1j*rng.normal(size=len(q))
    angle, shift = .18, .025
    result = higher_order_particle_remainders(k, q, c, angle, shift, domain=domain)
    assert result['degree_two_absolute_difference'] < 1e-13
    xyz, mass = cube_quadrature(20)
    for _ in range(3):
        u = rng.normal(size=5); u /= np.linalg.norm(u)
        a, b, d = angle*u[:3]
        omega = np.array([[0., -d, b], [d, 0., -a], [-b, a, 0.]])
        translation_phase = 2*np.pi*(q@(shift*u[3:]))
        nominal = np.exp(2j*np.pi*k@xyz.T)
        exact = np.real(c@np.exp(2j*np.pi*(k@expm(omega))@xyz.T+1j*translation_phase[:, None]))
        # Independently expand the scalar phase exponential using its Taylor
        # coefficient recurrence; no seminorm or bound routine is used here.
        phase_coefficients = [None]
        power = np.eye(3)
        for j in range(1, 6):
            power = power@omega
            derivative = 2*np.pi*(k@power)@xyz.T
            if j == 1:
                derivative += translation_phase[:, None]
            phase_coefficients.append(1j*derivative/math.factorial(j))
        exponential_coefficients = [np.ones_like(nominal)]
        approximation = np.real(c@nominal)
        for n in range(1, 6):
            coefficient = sum(j*phase_coefficients[j]*exponential_coefficients[n-j]
                              for j in range(1, n+1))/n
            exponential_coefficients.append(coefficient)
            approximation += np.real(c@(nominal*coefficient))
            if n >= 2:
                actual = np.sqrt(mass@((exact-approximation)**2))
                bound = next(r['field_remainder'] for r in result['records'] if r['taylor_degree'] == n)
                assert actual <= bound+3e-13


def test_invalid_degree_and_negative_bounds_are_rejected():
    with pytest.raises(ValueError):
        partial_bell([1., -2.])
    with pytest.raises(ValueError):
        higher_order_particle_remainders(None, None, None, .1, .01, degrees=[2, 2])
