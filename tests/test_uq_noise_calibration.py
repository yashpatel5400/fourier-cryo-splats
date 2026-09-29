import numpy as np
from scipy.stats import chi2
from fourier_splats.uq_noise_calibration import common_covariance_trace_upper


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
