import numpy as np
from fourier_splats.uq_projected_noise import fixed_noise_contrasts, projected_grouped_variance_upper
from fourier_splats.uq_noise_calibration import common_covariance_trace_upper


def test_fixed_family_projection_identities_and_rank_clamping():
    rng=np.random.default_rng(953101); x=rng.normal(size=(25,8)); pilot=rng.normal(size=(25,32))
    family,diag=fixed_noise_contrasts(x,pilot)
    assert list(family)==['raw','mean','ctf_4','ctf_8','pilot_4','pilot_8','pilot_16']
    assert [len(q) for q in family.values()]==[25,24,21,17,21,17,9]
    for name,q in family.items():
        np.testing.assert_allclose(q@q.T,np.eye(len(q)),atol=1e-13)
        if name!='raw': np.testing.assert_allclose(q@np.ones(25),0,atol=1e-13)
    rankone=np.ones((25,10))+np.arange(25)[:,None]
    _,small=fixed_noise_contrasts(rankone,rankone)
    assert small['ctf_8']['removed_rank']==2 and small['pilot_16']['removed_rank']==2


def test_simultaneous_selection_with_unequal_means_and_colored_noise():
    rng=np.random.default_rng(953102); n=25; dim=4; reps=7000
    design=rng.normal(size=(n,6));family,_=fixed_noise_contrasts(design,rng.normal(size=(n,18)))
    a=np.array([[1.,0.,0.,0.],[.7,.3,0.,0.],[.2,-.1,.2,0.],[-.3,.2,.1,.3]])
    sigma=a@a.T;v=rng.normal(size=(3,dim));groups=np.array([0,0,1]);sizes=np.array([2,2,1])
    inflated=v*np.sqrt(sizes[:,None]);truth=float(np.trace(inflated@sigma@inflated.T))
    means=design[:,:2]@rng.normal(size=(2,dim))*.7+3
    y=means+rng.normal(size=(reps,n,dim))@a.T
    expected=[]
    # Independent batched energy formula uses sample-space projectors and the
    # known covariance. It does not call the production selection function.
    weighted=y@inflated.T
    for q in family.values():
        t=common_covariance_trace_upper(np.zeros((len(q),1)),.05/7)['lower_tail_fraction']
        energy=np.einsum('bnk,nm,bmk->b',weighted,q.T@q,weighted)
        expected.append(energy/(len(q)*t))
    upper=np.min(expected,axis=0)
    assert np.mean(upper<truth)<.05+4*np.sqrt(.05*.95/reps)
    production=projected_grouped_variance_upper(y[0],v,groups,family,.05)
    np.testing.assert_allclose(production['estimator_variance_upper'],upper[0],rtol=1e-12)
    assert production['family_size']==7 and production['individual_failure_probability']==.05/7


def test_nonorthogonal_contrast_is_rejected():
    import pytest
    with pytest.raises(AssertionError):
        projected_grouped_variance_upper(np.zeros((4,2)),np.ones((1,2)),np.array([0]),{'bad':np.ones((2,4))})
