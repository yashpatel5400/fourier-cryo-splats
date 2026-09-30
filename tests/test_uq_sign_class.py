import numpy as np
import pytest
from fourier_splats.uq_sign_class import sign_class_support, one_sided_projection_support, projected_sign_bounds


@pytest.mark.parametrize('fraction', [0., .1, .5, 1.])
def test_support_matches_independent_conic_program(fraction):
    import cvxpy as cp
    h = np.array([1.7, -.5, .2, -2.1, .4]); r = 1.3; eta = fraction*r
    f = cp.Variable(len(h))
    program = cp.Problem(cp.Maximize(h@f), [cp.norm(f) <= r, cp.norm(cp.pos(-f)) <= eta])
    value = program.solve(solver='CLARABEL', tol_gap_abs=1e-10, tol_feas=1e-10, tol_gap_rel=1e-10)
    exact = sign_class_support(np.linalg.norm(np.maximum(h, 0)), np.linalg.norm(np.minimum(h, 0)), r, eta)
    np.testing.assert_allclose(exact, value, rtol=1e-8, atol=1e-8)


def test_continuous_linear_field_projection_converges_from_above():
    # h(x,y,z)=x on the unit cube: ||h||^2=1/12 and ||h_+||^2=1/24.
    truth = np.sqrt(1/24); bounds = []
    for m in [2, 4, 8, 16, 32]:
        x = (np.arange(m)+.5)/m-.5
        coefficients = np.broadcast_to(x[:, None, None], (m, m, m)).ravel()/m**1.5
        got = projected_sign_bounds(coefficients, np.sqrt(1/12), 1., 0.)
        assert got['upper_support'] >= truth-1e-14
        bounds.append(got['upper_support'])
    assert np.all(np.diff(bounds) < 0)
    assert bounds[-1]/truth-1 < .001


def test_projection_error_bound_covers_coefficient_perturbations():
    rng = np.random.default_rng(630934); h = rng.normal(size=13); error = rng.normal(size=13)*.03
    r, eta = 1.2, .3
    got = projected_sign_bounds(h+error, np.linalg.norm(h), r, eta, np.linalg.norm(error))
    actual = sign_class_support(np.linalg.norm(np.maximum(h, 0)), np.linalg.norm(np.minimum(h, 0)), r, eta)
    assert got['upper_support'] >= actual


def test_unconstrained_signs_recover_origin_ball():
    got = projected_sign_bounds(np.array([.2, -.3]), 1., 2., 2.)
    assert got['upper_support'] == got['negative_support_magnitude'] == got['bias_half_width'] == 2.
    assert got['center_offset'] == 0.


def test_invalid_projection_fails_closed():
    with pytest.raises(ValueError, match='Projection contradicts'):
        projected_sign_bounds(np.array([3., -4.]), 4., 1., .1)
    with pytest.raises(ValueError):
        one_sided_projection_support(1., 2., 1., .1)
