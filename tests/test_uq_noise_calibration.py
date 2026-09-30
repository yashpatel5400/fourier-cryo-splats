import numpy as np
from scipy.stats import chi2
from fourier_splats.uq_noise_calibration import common_covariance_trace_upper
from fourier_splats.uq_noise_calibration import estimator_variance_upper, uncenter_fourier_weights


def test_covariance_trace_chernoff_bound_against_exact_rank_one_tail():
    for n in [2,10,100,1000]:
        for beta in [.05,.005,1e-8]:
            result=common_covariance_trace_upper(np.ones((n,3)),beta)
            t=result['lower_tail_fraction']
            np.testing.assert_allclose(n*(t-1-np.log(t))/2,np.log(1/beta),rtol=1e-11)
            assert chi2.cdf(n*t,n)<=beta
            np.testing.assert_allclose(result['covariance_trace_upper'],3/t)


def test_covariance_trace_bound_is_invariant_to_rowwise_orthogonal_transforms():
    rng=np.random.default_rng(609671);y=rng.normal(size=(15,4));rotated=[]
    for row in y:
        q,_=np.linalg.qr(rng.normal(size=(4,4)));rotated.append(q@row)
    first=common_covariance_trace_upper(y);second=common_covariance_trace_upper(rotated)
    np.testing.assert_allclose(first['covariance_trace_upper'],second['covariance_trace_upper'],rtol=1e-14)


def test_estimator_variance_identity_with_correlated_fourier_coordinates():
    rng=np.random.default_rng(609852)
    n,m=7,10;v=rng.normal(size=(n,m));a=rng.normal(size=(m,m));sigma=a@a.T
    # Explicit stacked covariance verifies that projected calibration trace is
    # exactly the variance of the scalar independent-particle estimator.
    covariance=np.kron(np.eye(n),sigma)
    np.testing.assert_allclose(v.ravel()@covariance@v.ravel(),np.trace(v@sigma@v.T),rtol=1e-14)
    y=rng.normal(size=(31,m));res=estimator_variance_upper(y,v)
    expected=sum(float(row@v.T@v@row) for row in y)/(len(y)*res['calibration']['lower_tail_fraction'])
    np.testing.assert_allclose(res['estimator_variance_upper'],expected,rtol=1e-14)
    # With identical calibration data the projected energy cannot exceed the
    # old trace envelope times the squared Frobenius norm of V.
    old=common_covariance_trace_upper(y)['covariance_trace_upper']*np.sum(v*v)
    assert res['estimator_variance_upper']<=old


def test_centered_fourier_estimator_equals_raw_coordinate_pullback():
    rng=np.random.default_rng(609853);z=rng.normal(size=(11,9))+1j*rng.normal(size=(11,9))
    phases=-np.exp(1j*rng.normal(size=z.shape));w=rng.normal(size=(11,18));scale=2.7
    raw=uncenter_fourier_weights(w,phases,scale);centered=z*phases/scale
    lhs=np.sum(w*np.concatenate([centered.real,centered.imag],axis=1),axis=1)
    rhs=np.sum(raw*np.concatenate([z.real,z.imag],axis=1),axis=1)
    np.testing.assert_allclose(lhs,rhs,rtol=1e-14,atol=1e-14)


def test_projected_rank_one_calibration_failure_is_below_budget():
    # Exact distribution at arbitrary rank-one Sigma, correlated projected
    # coordinates and unequal nonzero means. Noncentrality only lowers failure.
    from scipy.stats import ncx2
    rng=np.random.default_rng(609854);n=23;v=rng.normal(size=(5,8));direction=rng.normal(size=8)
    signal=np.linspace(-.7,1.3,n);y=signal[:,None]*direction
    res=estimator_variance_upper(y,v,.005);t=res['calibration']['lower_tail_fraction']
    variance=float(np.sum((v@direction)**2))
    assert variance>0
    assert ncx2.cdf(n*t,n,float(signal@signal))<.005
