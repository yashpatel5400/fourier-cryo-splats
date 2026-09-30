import numpy as np
import pytest
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_shift_only import ShiftAwarePolynomialPoseFieldOperator


@pytest.mark.parametrize('backend', ['direct', 'nufft'])
@pytest.mark.parametrize('angle', [0., .01])
def test_zero_rotation_transform_reduction_preserves_field_adjoint_and_weight_derivative(backend, angle):
    rng = np.random.default_rng(3123); k = rng.normal(size=(3, 5, 3)); q = rng.normal(size=(3, 5, 2))
    ctf = rng.normal(size=(3, 5)); w = rng.normal(size=30)
    original = PolynomialPoseFieldOperator(k, q, ctf, w, 1.2, angle, .03, order=7, backend=backend, nthreads=1)
    optimized = ShiftAwarePolynomialPoseFieldOperator(k, q, ctf, w, 1.2, angle, .03, order=7, backend=backend, nthreads=1)
    scales = np.exp(rng.normal(size=60)); original.denominators = scales; optimized.denominators = scales
    u = rng.normal(size=60); v = rng.normal(size=343)
    np.testing.assert_allclose(optimized.matvec(u), original.matvec(u), rtol=1e-11, atol=1e-12)
    np.testing.assert_allclose(optimized.rmatvec(v), original.rmatvec(v), rtol=1e-11, atol=1e-12)
    np.testing.assert_allclose(optimized.weight_gradient(u, v), original.weight_gradient(u, v), rtol=1e-11, atol=1e-12)
    np.testing.assert_allclose(u@optimized.rmatvec(v), v@optimized.matvec(u), atol=1e-11)
