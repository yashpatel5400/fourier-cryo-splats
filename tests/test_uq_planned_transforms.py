import numpy as np
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_planned_transforms import PlannedQuadratureObservationGram,PlannedPolynomialPoseFieldOperator


def test_reused_gram_plans_match_reference_after_changing_weights():
    rng=np.random.default_rng(610321);k=rng.normal(size=(3,5,3));ctf=rng.normal(size=(3,5))
    reference=QuadratureObservationGram(k,ctf,.7,order=12,preconditioner_rank=0)
    planned=PlannedQuadratureObservationGram(k,ctf,.7,order=12,preconditioner_rank=0)
    for _ in range(3):
        w=rng.normal(size=30)
        np.testing.assert_allclose(planned.matvec(w),reference.matvec(w),rtol=1e-11,atol=1e-11)
    assert len(planned._transform_plans)==2


def test_reused_pose_plans_match_direct_fields_after_changing_weights():
    rng=np.random.default_rng(610322);k=rng.normal(size=(2,3,3));q=rng.normal(size=(2,3,2));ctf=rng.normal(size=(2,3))
    w=rng.normal(size=12)
    direct=PolynomialPoseFieldOperator(k,q,ctf,w,.8,.04,.003,order=10,backend='direct')
    planned=PlannedPolynomialPoseFieldOperator(k,q,ctf,w,.8,.04,.003,order=10)
    for _ in range(2):
        w=rng.normal(size=12);direct.set_weights(w);planned.set_weights(w)
        u=rng.normal(size=40);v=rng.normal(size=1000)
        np.testing.assert_allclose(planned.matvec(u),direct.matvec(u),rtol=1e-9,atol=1e-11)
        np.testing.assert_allclose(planned.rmatvec(v),direct.rmatvec(v),rtol=1e-9,atol=1e-11)
    assert len(planned._transform_plans)==2
