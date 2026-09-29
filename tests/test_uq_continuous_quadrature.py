import numpy as np
from fourier_splats.uq_continuous import ContinuousObservationGram,continuous_certificate
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram


def test_quadrature_gram_error_and_continuous_dual():
    rng=np.random.default_rng(609370);k=rng.normal(size=(2,5,3));ctf=rng.normal(size=(2,5));noise=.8;w=rng.normal(size=20)
    dense=ContinuousObservationGram(k,ctf,noise)
    # Coarse quadrature has a visible deterministic error; its uniform bound
    # must dominate both the scalar norm and observation-action discrepancies.
    low=QuadratureObservationGram(k,ctf,noise,order=3);error=low.quadrature_error(w);difference=low.matvec(w)-dense.matvec(w)
    assert abs(w@difference)<=error['squared_field_norm']+1e-8
    assert np.linalg.norm(difference)<=error['gram_action_norm']+1e-8
    high=QuadratureObservationGram(k,ctf,noise,order=32)
    np.testing.assert_allclose(high.matvec(w),dense.matvec(w),rtol=1e-9,atol=1e-9)
    left=continuous_certificate(dense,[[0,0,0]],[1],.16,.7,rtol=1e-5,maxiter=200)
    right=continuous_certificate(high,[[0,0,0]],[1],.16,.7,rtol=1e-5,maxiter=200)
    assert left['converged'] and right['converged']
    np.testing.assert_allclose(right['objective'],left['objective'],rtol=1e-7)
