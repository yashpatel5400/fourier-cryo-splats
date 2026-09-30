import numpy as np
from fourier_splats.uq_cubic_design import DifferentiableCubicPoseFieldOperator
from fourier_splats.uq_cubic_pose import gaussian_moments_cubic


def test_cubic_weight_gradient_directional_and_homogeneous():
    rng = np.random.default_rng(135)
    k = rng.normal(size=(2, 4, 3)); q = rng.normal(size=(2, 4, 2)); ctf = rng.normal(size=(2, 4))
    w = rng.normal(size=16)
    op = DifferentiableCubicPoseFieldOperator(k, q, ctf, w, 1.3, .12, .05, order=6, backend='direct')
    op.denominators = np.exp(rng.normal(size=110))
    u = rng.normal(size=110); v = rng.normal(size=op.shape[0]); dw = rng.normal(size=16)
    gradient = op.weight_gradient(u, v)
    np.testing.assert_allclose(w@gradient, v@op.matvec(u), atol=2e-13)
    step = 1e-5; op.set_weights(w+step*dw); plus = v@op.matvec(u)
    op.set_weights(w-step*dw); minus = v@op.matvec(u)
    np.testing.assert_allclose((plus-minus)/(2*step), dw@gradient, rtol=1e-9, atol=1e-10)


def test_pilot_moment_weight_gradient_matches_pairing_derivative():
    rng = np.random.default_rng(211)
    k = rng.normal(size=(2, 5, 3)); q = rng.normal(size=(2, 5, 2)); w = rng.normal(size=20)
    op = DifferentiableCubicPoseFieldOperator(k, q, np.ones((2, 5)), w, 2., .1, .03, order=5, backend='direct')
    moments = gaussian_moments_cubic(k, [[.1, .05, -.03]], [1.], .17)
    u = rng.normal(size=110); dw = rng.normal(size=20)
    gradient = op.moment_weight_gradient(u, moments)
    np.testing.assert_allclose(w@gradient, u@op.pair_moments(moments), atol=1e-13)
    step = 1e-5; op.set_weights(w+step*dw); plus = u@op.pair_moments(moments)
    op.set_weights(w-step*dw); minus = u@op.pair_moments(moments)
    np.testing.assert_allclose((plus-minus)/(2*step), dw@gradient, rtol=1e-9, atol=1e-10)
