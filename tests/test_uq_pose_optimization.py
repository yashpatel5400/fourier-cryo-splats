import numpy as np
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_pose_optimization import amplitude_gradient,cubic_penalty_coefficients,gram_spectral_power
from fourier_splats.uq_continuous_moments import integrated_cubic_remainder
from fourier_splats.uq_continuous_pose import polynomial_kernel_error


def test_shared_spectral_gradient_and_supporting_inequality():
    rng=np.random.default_rng(609621);n,nq=2,3
    k=.2*rng.normal(size=(n,nq,3));q=rng.normal(size=(n,nq,2))
    ctf=rng.normal(size=(n,nq));w=rng.normal(size=2*n*nq)
    op=PolynomialPoseFieldOperator(k,q,ctf,w,.7,.07,.002,order=4,backend='direct')
    op.establish_coefficient_scaling();coeff=op.coefficient_bound_matrix();error=polynomial_kernel_error(k,4)
    def evaluate(weights,with_gradient=False):
        op.set_weights(weights)
        f=np.column_stack([op.matvec(u) for u in np.eye(20*n)])
        left,singular,right=np.linalg.svd(f,full_matrices=False)
        wr=weights.reshape(n,2*nq);amp=np.hypot(wr[:,:nq],wr[:,nq:])
        masses=np.einsum('naj,nj->na',coeff,amp)
        value=np.sqrt(singular[0]**2+error*np.sum(masses**2))
        if not with_gradient:return value
        g=op.weight_gradient(right[0],left[:,0])
        gm=amplitude_gradient(np.einsum('naj,na->nj',coeff,masses),weights)
        return value,(singular[0]*g+error*gm)/value
    value,gradient=evaluate(w,True)
    direction=rng.normal(size=w.shape);step=1e-5
    finite=(evaluate(w+step*direction)-evaluate(w-step*direction))/(2*step)
    np.testing.assert_allclose(finite,gradient@direction,rtol=1e-7,atol=1e-9)
    np.testing.assert_allclose(gradient@w,value,rtol=1e-12)
    for _ in range(5):
        other=rng.normal(size=w.shape)
        assert evaluate(other)>=value+gradient@(other-w)-1e-10


def test_particle_specific_cubic_penalty_matches_original_scalar_formula():
    rng=np.random.default_rng(609622);n,nq=3,5
    k=rng.normal(size=(n,nq,3));q=rng.normal(size=(n,nq,2));ctf=rng.normal(size=(n,nq))
    w=rng.normal(size=(n,2*nq));angles=np.array([.01,.04,.08]);shifts=np.array([.001,.002,.003])
    coef=cubic_penalty_coefficients(k,q,ctf/.8,angles,shifts,3.)
    new=np.sum(coef*np.hypot(w[:,:nq],w[:,nq:]))
    old=sum(integrated_cubic_remainder(k[i:i+1],q[i:i+1],ctf[i:i+1],w[i:i+1],.8,
                                     angles[i],shifts[i],2.,1.)['frequency_weighted_remainder_bias'] for i in range(n))
    np.testing.assert_allclose(new,old,rtol=1e-12)


def test_low_rank_gram_power_against_dense_symmetric_power():
    from types import SimpleNamespace
    rng=np.random.default_rng(609623);n,nq=2,3;factors=[];matrices=[]
    for _ in range(2):
        q,_=np.linalg.qr(rng.normal(size=(n*nq,3)))
        lam=np.array([2.,5.,11.]);f=q*np.sqrt(lam)
        factors.append((f,lam));matrices.append(f@f.T+3*np.eye(n*nq))
    pack=lambda r,i:np.concatenate([r.reshape(n,nq),i.reshape(n,nq)],axis=1).ravel()
    gram=SimpleNamespace(n=n,nq=nq,preconditioner_factors=factors,pack=pack)
    vector=rng.normal(size=2*n*nq);v=vector.reshape(n,2*nq)
    for power in [-.5,.5]:
        expected=[]
        for matrix,part in zip(matrices,[v[:,:nq].ravel(),v[:,nq:].ravel()]):
            e,u=np.linalg.eigh(matrix);expected.append((u*e**power)@(u.T@part))
        np.testing.assert_allclose(gram_spectral_power(gram,vector,3.,power),pack(*expected),rtol=1e-12)
