import numpy as np
from scipy.linalg import block_diag
from types import SimpleNamespace
from fourier_splats.uq_block_preconditioner import BlockLowRankCoordinates, cubic_remainder_metric


def test_dense_whitening_inverse_and_nonsymmetric_gradient_map():
    rng = np.random.default_rng(650151); n, q, rank = 3, 4, 5
    a = rng.normal(size=(n, 2, q, q)); blocks = a@a.swapaxes(-1, -2)+.2*np.eye(q)
    factors = [rng.normal(size=(n*q, rank)) for _ in range(2)]
    metric = block_diag(*blocks.reshape(n*2, q, q))
    for coordinate, factor in enumerate(factors):
        indices = np.array([i*2*q+coordinate*q+j for i in range(n) for j in range(q)])
        metric[np.ix_(indices, indices)] += factor@factor.T
    transform = BlockLowRankCoordinates(blocks, factors)
    eye = np.eye(2*n*q)
    T = np.column_stack([transform.to_weights(e) for e in eye])
    np.testing.assert_allclose(T.T@metric@T, eye, atol=1e-12)
    np.testing.assert_allclose(np.column_stack([transform.to_coordinates(e) for e in eye])@T, eye, atol=1e-12)
    assert np.linalg.norm(T-T.T) > .01
    gradient = rng.normal(size=len(eye))
    np.testing.assert_allclose(transform.gradient_to_coordinates(gradient), T.T@gradient, atol=1e-12)
    x = rng.normal(size=len(eye)); direction = rng.normal(size=x.shape)
    f = lambda z: np.logaddexp(0, transform.to_weights(z)).sum()
    w = transform.to_weights(x); g = 1/(1+np.exp(-w))
    step = 1e-5
    np.testing.assert_allclose((f(x+step*direction)-f(x-step*direction))/(2*step),
        transform.gradient_to_coordinates(g)@direction, rtol=1e-8, atol=1e-9)


def test_local_remainder_metric_matches_its_independent_dense_hessian_majorant():
    rng = np.random.default_rng(650152); n, q, d = 2, 3, 4
    a = rng.normal(size=(n, d, 2, q, q)); blocks = a@a.swapaxes(-1, -2)+.1*np.eye(q)
    multipliers = rng.uniform(.1, 1., size=(n, d)); w = rng.normal(size=(n, 2, q))
    factors = [rng.normal(size=(n*q, 2)) for _ in range(2)]
    gram = SimpleNamespace(preconditioner_factors=[(f, np.linalg.eigvalsh(f.T@f)) for f in factors])
    penalty = SimpleNamespace(n=n, nq=q, blocks=blocks, multipliers=multipliers)
    transform = cubic_remainder_metric(gram, penalty, w.ravel(), 2., 1., .3, 2.7)
    expected = 2.7/np.linalg.norm(w)*np.eye(2*n*q)
    for i in range(n):
        sl = slice(i*2*q, (i+1)*2*q)
        for j in range(d):
            G = block_diag(*blocks[i, j]); value = np.sqrt(w[i].ravel()@G@w[i].ravel())
            expected[sl, sl] += 3*multipliers[i, j]/value*G
    for c, f in enumerate(factors):
        indices = np.array([i*2*q+c*q+j for i in range(n) for j in range(q)])
        expected[np.ix_(indices, indices)] += (2/.3)*f@f.T
    T = np.column_stack([transform.to_weights(e) for e in np.eye(2*n*q)])
    np.testing.assert_allclose(T.T@expected@T, np.eye(2*n*q), atol=1e-12)


def test_zero_rank_and_degenerate_low_rank_factors():
    blocks = np.broadcast_to(np.diag([1., 4.]), (2, 2, 2, 2)).copy()
    for factors in [[np.empty((4, 0))]*2, [np.zeros((4, 3))]*2]:
        transform = BlockLowRankCoordinates(blocks, factors)
        w = np.arange(8, dtype=float)
        np.testing.assert_allclose(transform.to_weights(transform.to_coordinates(w)), w, atol=1e-14)


def test_zero_derivative_grams_have_finite_noise_only_metric():
    penalty = SimpleNamespace(n=2, nq=3, blocks=np.zeros((2, 4, 2, 3, 3)), multipliers=np.ones((2, 4)))
    gram = SimpleNamespace(preconditioner_factors=[])
    w = np.arange(12, dtype=float)
    metric = cubic_remainder_metric(gram, penalty, w, 2., 1., .3, 2.7)
    np.testing.assert_allclose(metric.to_weights(w), w/np.sqrt(2.7/np.linalg.norm(w)), atol=1e-13)
    assert metric.diagnostics['retained_ranks'] == [0, 0]
