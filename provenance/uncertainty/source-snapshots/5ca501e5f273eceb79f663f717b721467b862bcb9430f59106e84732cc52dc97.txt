import numpy as np
from fourier_splats.uq_cached_quadrature import CachedQuadratureGram
from fourier_splats.uq_continuous import ContinuousObservationGram
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram


def test_cached_transform_matches_uncached_and_exact_continuous_kernel():
    rng = np.random.default_rng(927)
    k = rng.uniform(-2, 2, (3, 9, 3)); ctf = rng.normal(size=(3, 9))
    exact = ContinuousObservationGram(k, ctf, .5)
    base = QuadratureObservationGram(k, ctf, .5, order=24, preconditioner_rank=0)
    cached = CachedQuadratureGram(base)
    for _ in range(3):
        weights = rng.normal(size=54)
        np.testing.assert_allclose(cached.field(weights), base.field(weights), rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(cached.matvec(weights), exact.matvec(weights), rtol=1e-9, atol=1e-9)
        np.testing.assert_allclose(cached.matvec(weights), base.matvec(weights), rtol=1e-10, atol=1e-10)
