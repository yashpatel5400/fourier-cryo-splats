import numpy as np
from types import SimpleNamespace
from fourier_splats.uq_noise_design import CenteredNoiseMetric,NoiseMetricGram,shrunk_second_moment


def test_noise_metric_powers_and_colored_variance_against_explicit_blocks():
    rng=np.random.default_rng(609871);n,nq=3,4;phases=np.exp(1j*rng.normal(size=(n,nq)))
    a=rng.normal(size=(2*nq,2*nq));cov=a@a.T+np.eye(2*nq)
    metric=CenteredNoiseMetric(cov,phases);w=rng.normal(size=n*2*nq)
    blocks=[]
    for p in phases:
        o=np.block([[np.diag(p.real),-np.diag(p.imag)],[np.diag(p.imag),np.diag(p.real)]])
        blocks.append(o@cov@o.T)
    from scipy.linalg import block_diag
    c=block_diag(*blocks);values,vectors=np.linalg.eigh(c)
    expected=(vectors*values**-.5)@vectors.T
    np.testing.assert_allclose(metric.apply(np.eye(len(w))),expected,rtol=1e-11,atol=1e-12)
    np.testing.assert_allclose(metric.apply(metric.apply(w),.5),w,rtol=1e-11,atol=1e-12)
    transformed=metric.apply(w);np.testing.assert_allclose(transformed@c@transformed,w@w,rtol=1e-12)


def test_metric_gram_and_woodbury_preconditioner_against_dense_linear_algebra():
    rng=np.random.default_rng(609872);n,nq=2,3;dim=2*n*nq
    factors=[];full=[]
    for block in range(2):
        f=rng.normal(size=(n*nq,3));factors.append((f,np.zeros(3)))
        packed=np.zeros((n,2*nq,3));packed[:,block*nq:(block+1)*nq]=f.reshape(n,nq,3);full.append(packed.reshape(dim,3))
    f=np.concatenate(full,axis=1);g=f@f.T;a=rng.normal(size=dim)
    gram=SimpleNamespace(shape=(dim,dim),n=n,nq=nq,preconditioner_factors=factors,
        matvec=lambda x:g@x,target=lambda *args:(a,2.),
        quadrature_error=lambda x:{'squared_field_norm':.1*np.linalg.norm(x)**2,'gram_action_norm':.2*np.linalg.norm(x)})
    metric=CenteredNoiseMetric(shrunk_second_moment(rng.normal(size=(17,2*nq))),np.exp(1j*rng.normal(size=(n,nq))))
    wrapped=NoiseMetricGram(gram,metric);w=metric.apply(np.eye(dim));x=rng.normal(size=dim)
    np.testing.assert_allclose(wrapped.matvec(x),w@g@w@x,rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(wrapped.target(None,None,None)[0],w@a,rtol=1e-12)
    np.testing.assert_allclose(wrapped.preconditioner(1.7)@x,np.linalg.solve(w@g@w+1.7*np.eye(dim),x),rtol=1e-11,atol=1e-12)
    errors=wrapped.quadrature_error(x);base=gram.quadrature_error(w@x)
    np.testing.assert_allclose(errors['squared_field_norm'],base['squared_field_norm'],rtol=1e-12)
    np.testing.assert_allclose(errors['gram_action_norm'],np.linalg.norm(w,2)*base['gram_action_norm'],rtol=1e-12)
