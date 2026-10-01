"""Independent block-Gaussian checks, not a dataset coverage experiment."""
import numpy as np
import pytest
from scipy.linalg import block_diag
from fourier_splats.uq_paired_power import (complex_log_normalizer,
    common_mean_coefficient, mismatch_coefficient, paired_power_log_factor)


def realify(z):
    return np.column_stack([z.real,z.imag]).ravel()


def direct_log_mgf(mean_a,mean_b,cov_a,cov_b,weights):
    n=len(mean_a);t=np.diag(np.repeat(weights,2));zero=np.zeros_like(t)
    h=np.block([[zero,t/2],[t/2,zero]])
    cov=block_diag(cov_a,cov_b);m=np.r_[realify(mean_a),realify(mean_b)]
    values,vectors=np.linalg.eigh(cov)
    root=(vectors*np.sqrt(np.maximum(values,0)))@vectors.T
    precision=np.eye(4*n)-2*root@h@root
    assert np.linalg.eigvalsh(precision).min()>0
    linear=2*root@h@m
    return -.5*np.linalg.slogdet(precision)[1]+m@h@m+.5*linear@np.linalg.solve(precision,linear)


def test_complex_normalizer_and_noncentral_factor_match_dense_gaussian():
    rng=np.random.default_rng(261005)
    for _ in range(25):
        v=rng.uniform(.2,3.,4);t=rng.uniform(-.8,.8,4)/v
        mu=rng.normal(size=4)+1j*rng.normal(size=4)
        d=np.diag(np.repeat(v,2))
        exact=direct_log_mgf(mu,mu,d,d,t)
        expected=-complex_log_normalizer(t,v).sum()+np.dot(common_mean_coefficient(t,v),abs(mu)**2)
        np.testing.assert_allclose(exact,expected,rtol=2e-12,atol=2e-12)


def test_covariance_upper_bound_with_noncommuting_correlations_and_mismatch():
    rng=np.random.default_rng(261006)
    for _ in range(100):
        n=4;v=rng.uniform(.2,3.,n);t=rng.uniform(-.85,.85,n)/v
        droot=np.diag(np.sqrt(np.repeat(v,2)))
        covs=[]
        for _ in range(2):
            a=rng.normal(size=(2*n,2*n));s=a@a.T
            s*=rng.uniform(.1,.99)/np.linalg.eigvalsh(s).max()
            covs.append(droot@s@droot)
        mu=rng.normal(size=n)+1j*rng.normal(size=n)
        lo=rng.uniform(.5,.95,n);hi=rng.uniform(1.05,1.5,n);rad=rng.uniform(0,np.pi,n)
        gain=rng.uniform(lo,hi);phase=rng.uniform(-rad,rad)
        second=mu*gain*np.exp(1j*phase)
        direct=direct_log_mgf(mu,second,*covs,t)
        bound=-complex_log_normalizer(t,v).sum()+np.dot(mismatch_coefficient(t,v,lo,hi,rad),abs(mu)**2)
        assert direct<=bound+2e-11


def test_mismatch_endpoints_attain_bound_at_bounding_covariance():
    t=np.array([.3,-.25,.1,-.15]);v=np.array([1.,1.4,2.,.7]);lo=np.array([.5,.8,.9,.6]);hi=np.array([1.2,1.1,1.5,1.3]);rad=np.array([.3,1.,.4,4.])
    mu=np.array([1+.2j,2-1j,-1+.1j,.7+1j]);d=np.diag(np.repeat(v,2))
    phase=np.where(t>=0,0.,np.minimum(rad,np.pi))
    c_lo=mismatch_coefficient(t,v,lo,lo,rad);c_hi=mismatch_coefficient(t,v,hi,hi,rad)
    gain=np.where(c_lo>c_hi,lo,hi)
    expected=-complex_log_normalizer(t,v).sum()+np.dot(mismatch_coefficient(t,v,lo,hi,rad),abs(mu)**2)
    np.testing.assert_allclose(direct_log_mgf(mu,mu*gain*np.exp(1j*phase),d,d,t),expected,rtol=2e-12,atol=2e-12)


def test_zero_variance_and_identical_mean_reduction():
    t=np.array([.1,-.2]);v=np.array([2.,3.])
    np.testing.assert_allclose(mismatch_coefficient(t,v),common_mean_coefficient(t,v))
    np.testing.assert_allclose(mismatch_coefficient(t,0),t)
    np.testing.assert_equal(complex_log_normalizer(t,0),np.zeros(2))


def test_factor_uses_complex_real_coordinate_variance_convention():
    x=np.array([1+2j,-1+.5j]);z=np.array([.5-.2j,3+1j]);t=np.array([.2,-.1]);v=np.array([2.,3.])
    direct=np.sum(t*(x.real*z.real+x.imag*z.imag)+np.log(1-(t*v)**2))-.7
    assert paired_power_log_factor(x,z,t,v,.7)==pytest.approx(direct)
    with pytest.raises(ValueError):paired_power_log_factor(x,z,t,v,-.1)
    with pytest.raises(ValueError):complex_log_normalizer([1.],[1.])
    with pytest.raises(ValueError):mismatch_coefficient(t,v,gain_low=2,gain_high=1)
    with pytest.raises(ValueError):mismatch_coefficient(t,v,phase_radius=-1)


def test_unbounded_phase_mismatch_removes_signed_signal_cancellation():
    t=np.array([.1,-.1]);v=np.ones(2)
    assert np.all(mismatch_coefficient(t,v,phase_radius=np.pi)>0)
