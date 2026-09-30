import numpy as np
import pytest
from fourier_splats.uq_sobolev_penalty import derivative_gram_blocks, SobolevRemainderPenalty
from fourier_splats.uq_ball_remainder import ball_sobolev_norms
from fourier_splats.uq_higher_remainder import higher_order_particle_remainders


@pytest.mark.parametrize('domain', ['cube', 'ball'])
def test_derivative_grams_match_independent_complex_pair_identity(domain):
    rng = np.random.default_rng(913); k = rng.normal(size=(8, 3)); t = rng.normal(size=8)
    w = rng.normal(size=(2, 8)); orders = (0, 1, 2, 3, 4)
    grams, pads = derivative_gram_blocks(k, t, .6, orders=orders, domain=domain)
    values = np.sqrt(np.einsum('dcij,ci,cj->d', grams, w, w))
    expected, _ = ball_sobolev_norms(k, t*(w[0]+1j*w[1]), .6, orders=orders, domain=domain)
    np.testing.assert_allclose(values, expected, rtol=1e-10)
    assert all(p['diagonal_pad'] >= 0 for p in pads)
    assert np.min(np.linalg.eigvalsh(grams)) >= 0


@pytest.mark.parametrize('nonplanar', [False, True])
def test_remainder_subgradient_euler_convexity_and_independent_value(nonplanar):
    rng = np.random.default_rng(192); q = rng.normal(size=(2, 7, 2)); k = np.concatenate([q, np.zeros((2, 7, 1))], axis=-1)
    if nonplanar: k += rng.normal(size=k.shape)*.01
    transfer = rng.normal(size=(2, 7)); w = rng.normal(size=28)
    penalty = SobolevRemainderPenalty(k, q, transfer, .11, .03)
    value, gradient = penalty.value_gradient(w)
    expected = 0.
    for i in range(2):
        row = w.reshape(2, 2, 7)[i]
        c = transfer[i]*(row[0]+1j*row[1])
        expected += higher_order_particle_remainders(k[i], q[i], c, .11, .03, degrees=(3,), domain='cube')['records'][0]['field_remainder']
    np.testing.assert_allclose(value, expected, rtol=1e-9)
    np.testing.assert_allclose(w@gradient, value, rtol=1e-13)
    for _ in range(4):
        direction = rng.normal(size=w.shape); step = 1e-5
        finite = (penalty.value_gradient(w+step*direction)[0]-penalty.value_gradient(w-step*direction)[0])/(2*step)
        np.testing.assert_allclose(finite, gradient@direction, rtol=1e-8, atol=1e-10)
        other = w+direction
        assert penalty.value_gradient(other)[0] >= value+gradient@direction-1e-12
    zero, subgradient = penalty.value_gradient(np.zeros_like(w))
    assert zero == 0 and np.all(subgradient == 0)


def test_zero_transfer_and_invalid_inputs():
    q = np.array([[[1., 0.], [0., 1.], [1., 1.]]]); k = np.pad(q, ((0, 0), (0, 0), (0, 1)))
    penalty = SobolevRemainderPenalty(k, q, np.zeros((1, 3)), .1, .01)
    v, g = penalty.value_gradient(np.ones(6)); assert v == 0 and np.all(g == 0)
    with pytest.raises(ValueError): SobolevRemainderPenalty(k, q, np.ones((1, 3)), .1, .01, degree=0)
    with pytest.raises(ValueError): derivative_gram_blocks(k[0], np.ones(2), .5)
    with pytest.raises(ValueError): derivative_gram_blocks(k[0], np.ones(3), .5, domain='invalid')
