import numpy as np
from fourier_splats.physics import ctf
from fourier_splats.uq_ctf_uncertainty import ctf_uncertainty_envelope,continuous_ctf_bias_bound


def test_ctf_envelope_covers_simultaneous_nonlinear_parameter_changes():
    rng=np.random.default_rng(609691);q=rng.uniform(-.15,.15,(30,2))
    p=np.tile([360,1.34,14000,18000,17,300,2.7,.1,5.],(3,1))
    df=np.array([0,100,500]);angle=np.array([0,2,5]);phase=np.array([0,1,4])
    gain=np.array([0,.03,.1]);bf=np.array([0,10,30])
    audit=ctf_uncertainty_envelope(q,p,df,angle,phase,gain,bf)
    np.testing.assert_allclose(audit['nominal_transfer_float64'],ctf(q,p),atol=4e-8,rtol=0)
    np.testing.assert_allclose(audit['absolute_transfer_error'][0],0,atol=1e-15)
    for _ in range(100):
        other=p.copy();other[:,2:4]+=rng.uniform(-1,1,(3,2))*df[:,None]
        other[:,4]+=rng.uniform(-1,1,3)*angle;other[:,8]+=rng.uniform(-1,1,3)*phase
        factor=(1+rng.uniform(-1,1,(3,1))*gain[:,None])*np.exp(-rng.uniform(-1,1,(3,1))*bf[:,None]*np.sum(q*q,axis=1)[None]/4)
        actual=ctf(q,other)*factor
        assert np.all(actual>=audit['minimum_transfer']-5e-8)
        assert np.all(actual<=audit['maximum_transfer']+5e-8)


def test_ctf_bias_bound_covers_independent_direct_spatial_integral():
    rng=np.random.default_rng(609692);n,nq=3,5
    w=rng.normal(size=(n,2*nq));envelope=rng.uniform(.01,.2,(n,nq));noise=.7;R=2.
    xyz=rng.uniform(-.5,.5,(5000,3));k=rng.normal(size=(n,nq,3));shift=rng.normal(size=(n,nq))
    error=rng.uniform(-1,1,(n,nq))*envelope
    c=(w[:,:nq]+1j*w[:,nq:])*error/noise
    field=np.real(np.einsum('nq,nqx->x',c,np.exp(2j*np.pi*(np.einsum('nqc,xc->nqx',k,xyz)+shift[:,:,None]))))
    sampled_norm=np.sqrt(np.mean(field**2))
    assert R*sampled_norm<=continuous_ctf_bias_bound(w,noise,envelope,R)
