import numpy as np
from fourier_splats.uq_physics import (
    evaluate_pairs, image_design, pose_jacobian, perturbed_images,
    pose_remainder_bounds, pose_design_derivatives, realify,gaussian_pair_gram,density_energy_coordinates,
)


def example(seed=51):
    rng=np.random.default_rng(seed)
    k=rng.normal(size=(7,19,3))*2
    q=rng.normal(size=(7,19,2))*2
    centers=rng.normal(size=(11,3))*2
    sigma=rng.uniform(0.5,1.3,11)
    c=rng.normal(size=22)
    ctf=rng.uniform(-1,1,(7,19))
    return rng,k,q,centers,sigma,c,ctf


def test_realified_gaussian_forward_and_hermitian():
    _,k,_,centers,sigma,c,ctf=example()
    f=evaluate_pairs(k,centers,sigma,c)
    fm=evaluate_pairs(-k,centers,sigma,c)
    assert np.allclose(f.conjugate(),fm)
    a=image_design(k,ctf,centers,sigma)
    assert np.allclose(a@c,realify(f*ctf))


def test_pose_jacobian_against_finite_difference():
    rng,k,q,centers,sigma,c,ctf=example()
    jac=pose_jacobian(k,q,ctf,centers,sigma,c,32,0.05,0.3)
    for axis in range(5):
        u=np.zeros((len(k),5));u[:,axis]=1e-5
        fp=perturbed_images(k,q,ctf,centers,sigma,c,32,u,0.05,0.3)
        fm=perturbed_images(k,q,ctf,centers,sigma,c,32,-u,0.05,0.3)
        numeric=realify((fp-fm)/(2e-5))
        assert np.allclose(jac[:,:,axis],numeric,rtol=1e-6,atol=1e-8)


def test_uniform_remainder_bound_under_joint_density_pose_errors():
    rng,k,q,centers,sigma,c,ctf=example()
    B=0.4
    for angle in [0.001,0.05,0.3]:
        jac=pose_jacobian(k,q,ctf,centers,sigma,c,32,angle,0.4)
        a=image_design(k,ctf,centers,sigma)
        gamma,L,H=pose_remainder_bounds(k,q,ctf,centers,sigma,c,32,angle,0.4,B)
        for _ in range(40):
            delta=rng.normal(size=len(c));delta*=B/np.linalg.norm(delta)
            u=rng.normal(size=(len(k),5));u/=np.linalg.norm(u,axis=1,keepdims=True)
            truth=realify(perturbed_images(k,q,ctf,centers,sigma,c+delta,32,u,angle,0.4))
            remainder=truth-a@(c+delta)-np.einsum('nmq,nq->nm',jac,u)
            assert np.all(np.linalg.norm(remainder,axis=1)<=gamma+1e-10)


def test_structured_density_pose_interaction_and_remainder():
    rng,k,q,centers,sigma,c,ctf=example(71)
    B=.4;angle=.07;shift=.3;box=32
    derivatives=pose_design_derivatives(k,q,ctf,centers,sigma,box,angle,shift)
    jac=pose_jacobian(k,q,ctf,centers,sigma,c,box,angle,shift)
    assert np.allclose(np.einsum('nmpq,p->nmq',derivatives,c),jac)
    gamma,_,_=pose_remainder_bounds(k,q,ctf,centers,sigma,c,box,angle,shift,B,structured=True)
    a=image_design(k,ctf,centers,sigma)
    for _ in range(80):
        delta=rng.normal(size=len(c));delta*=B/np.linalg.norm(delta)
        u=rng.normal(size=(len(k),5));u/=np.linalg.norm(u,axis=1,keepdims=True)
        truth=realify(perturbed_images(k,q,ctf,centers,sigma,c+delta,box,u,angle,shift))
        linear=np.einsum('nmpq,p,nq->nm',derivatives,c+delta,u)
        remainder=truth-a@(c+delta)-linear
        assert np.all(np.linalg.norm(remainder,axis=1)<=gamma+1e-10)
        w=rng.normal(size=a.shape[:2])
        cross=np.einsum('nmpq,p,nq,nm->n',derivatives,delta,u,w)
        bound=B*np.linalg.norm(np.einsum('nmpq,nm->npq',derivatives,w),axis=(1,2))
        assert np.all(np.abs(cross)<=bound+1e-10)


def test_analytic_density_energy_against_independent_quadrature():
    rng=np.random.default_rng(415)
    centers=np.array([[0.,0.,0.],[1.2,.4,-.6],[-.8,.3,.7]])
    sigma=np.array([.5,.8,.65]);c=rng.normal(size=6)
    gram=gaussian_pair_gram(centers,sigma)
    g=np.arange(-6,6.01,.15);z,y,x=np.meshgrid(g,g,g,indexing='ij')
    k=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    numerical=np.sum(np.abs(evaluate_pairs(k,centers,sigma,c))**2)*.15**3
    assert np.isclose(c@gram@c,numerical,rtol=1e-8)
    transform,record=density_energy_coordinates(centers,sigma)
    assert record['excluded_modes']==1
    assert np.allclose(transform.T@gram@transform,np.eye(5),atol=1e-12)
