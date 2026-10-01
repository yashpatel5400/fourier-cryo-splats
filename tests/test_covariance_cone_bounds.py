import numpy as np
from scipy.optimize import minimize_scalar
from fourier_splats.uq_covariance_cone_bounds import (unrestricted_matrix_growth_upper,
    scalar_dual_tangent_upper,full_covariance_cone,log_growth_variance)


def test_spectral_upper_matches_independent_scalar_maxima():
    values=np.array([-4.,-.2,0.,.001,2.])
    actual=unrestricted_matrix_growth_upper(np.diag(values))
    reference=0.
    for d in values:
        fit=minimize_scalar(lambda t:-(d*t+.5*np.log1p(-t*t)),bounds=(-1+1e-10,1-1e-10),method='bounded',options={'xatol':1e-14})
        reference-=fit.fun
    np.testing.assert_allclose(actual['raw_spectral_value'],reference,rtol=1e-11,atol=1e-11)
    assert actual['upper']>=reference


def test_scalar_tangent_dominates_dense_grid_including_near_boundary():
    for s,a in [(1.,1.),(1e5,0.),(.01,2.),(100.,.001),(1e-12,1e-12)]:
        result=scalar_dual_tangent_upper(s,a)
        grid=np.concatenate([np.linspace(-.999,.999,2001),[result['root']],1-np.geomspace(1e-12,.001,100)])
        values=s*grid+np.log1p(-grid*grid)-a*grid/(1-grid)
        assert result['upper']>=values.max()-1e-9


def test_full_cone_exact_and_separated_cases():
    null=np.array([[1.,1.],[-1.,-1.]])
    exact=full_covariance_cone(null,np.ones((2,2)),maximum_seconds=10)
    assert exact['residual_frobenius']<1e-9
    assert exact['growth']['1.0']['upper']<1e-9
    different=full_covariance_cone(null,np.array([[1.,-1.],[-1.,1.]]),maximum_seconds=10)
    assert different['growth']['1.0']['expected_log_lower']>.1
    assert different['growth']['1.0']['gap']>=-1e-10


def test_log_variance_matches_direct_gaussian_simulation():
    rng=np.random.default_rng(651);m=np.array([[1.,2.],[-.2,.5],[0.,-1.]])
    t=np.array([[.2,.1],[.1,-.15]])
    theory=log_growth_variance(m,t)
    indices=rng.integers(0,3,200000);x=m[indices]+rng.normal(size=(len(indices),2));z=m[indices]+rng.normal(size=x.shape)
    statistic=np.einsum('ni,ij,nj->n',x,t,z)+.5*np.log1p(-np.linalg.eigvalsh(t)**2).sum()
    assert abs(statistic.mean()-theory['mean'])<.008
    assert abs(statistic.var()-theory['variance'])<.015
