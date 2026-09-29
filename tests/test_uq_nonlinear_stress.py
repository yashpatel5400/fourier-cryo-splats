import numpy as np
import torch
from fourier_splats.uq_data import VoxelObservationOperator,SupportedVoxelOperator
from fourier_splats.uq_pose_audit import pose_derivative_adjoints
from fourier_splats.uq_nonlinear_stress import TorchPoseAdjoint,exact_pose_adjoint,eliminated_bias


def test_nonlinear_adjoint_gradient_and_density_elimination():
    rng=np.random.default_rng(609351);n,q,box=3,5,6
    k=rng.normal(size=(n,q,3));detector=rng.normal(size=(n,q,2));ctf=rng.normal(size=(n,q))
    op=SupportedVoxelOperator(VoxelObservationOperator(k,ctf,box,.7),rng.uniform(size=box**3)>.4)
    w=rng.normal(size=(n,2*q));model=TorchPoseAdjoint(op,detector,w,.19,.23)
    u=torch.tensor(rng.normal(size=(n,5))*.2,dtype=torch.float64,requires_grad=True)
    np.testing.assert_allclose(model(u).detach(),exact_pose_adjoint(op,detector,w,u.detach(),.19,.23),atol=2e-13)
    assert torch.autograd.gradcheck(model,(u,),eps=1e-6,atol=1e-6,rtol=1e-5)
    zero=torch.zeros_like(u,requires_grad=True);pilot=torch.tensor(rng.normal(size=op.shape[1]))
    nominal=model(zero);np.testing.assert_allclose(nominal.detach(),op.adjoint(w.ravel()),atol=1e-8)
    grad=torch.autograd.grad(pilot@nominal,zero)[0].detach().numpy()
    for i in range(n):
        d1,_=pose_derivative_adjoints(op,detector,w,.19,.23,i)
        np.testing.assert_allclose(grad[i],pilot.numpy()@d1,atol=1e-10)
    ell=torch.tensor(rng.normal(size=op.shape[1]));B=.9
    a=model(u);h=a-ell
    for sign in [-1,1]:
        delta=sign*B*h/torch.linalg.vector_norm(h)
        actual=sign*(torch.dot(pilot,a-nominal)+torch.dot(delta,h))
        np.testing.assert_allclose(actual.detach(),eliminated_bias(a,nominal,ell,pilot,B,sign).detach(),atol=1e-12)
