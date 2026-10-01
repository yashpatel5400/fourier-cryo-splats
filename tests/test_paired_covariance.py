import numpy as np
import pytest
from fourier_splats.uq_paired_covariance import (matrix_gaussian_terms,
    matrix_log_factor, covariance_cone_diagnostic, direction_growth, svec, smat)


def dense_mgf(t, mean, first_cov, second_cov):
    n = len(t)
    h = np.block([[np.zeros_like(t),t/2],[t/2,np.zeros_like(t)]])
    cov = np.block([[first_cov,np.zeros_like(t)],[np.zeros_like(t),second_cov]])
    mu = np.tile(mean,2); operator = np.eye(2*n)-2*cov@h
    eigenvalues,vectors=np.linalg.eigh(cov)
    root=(vectors*np.sqrt(np.maximum(eigenvalues,0)))@vectors.T
    assert np.linalg.eigvalsh(np.eye(2*n)-2*root@h@root).min()>0
    sign, logdet = np.linalg.slogdet(operator)
    assert sign > 0
    return -.5*logdet+mu@h@np.linalg.solve(operator,mu)


def test_matrix_formula_and_covariance_monotonicity_against_dense_mgf():
    rng = np.random.default_rng(261004)
    for _ in range(100):
        q,_ = np.linalg.qr(rng.normal(size=(5,5)))
        t = (q*rng.uniform(-.8,.8,5))@q.T; mu = rng.normal(size=5)
        lognorm,u = matrix_gaussian_terms(t)
        exact = dense_mgf(t,mu,np.eye(5),np.eye(5))
        np.testing.assert_allclose(exact, -lognorm+mu@u@mu, rtol=1e-12,atol=1e-12)
        covariances=[]
        for _ in range(2):
            a=rng.normal(size=(5,5));c=a@a.T
            covariances.append(c/np.linalg.norm(c,2)*rng.uniform(.01,1))
        assert dense_mgf(t,mu,*covariances) <= exact+1e-11


def test_svec_preserves_inner_product_and_inversion():
    rng=np.random.default_rng(180)
    a=rng.normal(size=(7,7));a=(a+a.T)/2
    b=rng.normal(size=(7,7));b=(b+b.T)/2
    np.testing.assert_allclose(svec(a)@svec(b),np.sum(a*b))
    np.testing.assert_allclose(smat(svec(a),7),a)


def test_known_mixture_is_replayed_and_cannot_grow():
    rng=np.random.default_rng(912);m=rng.normal(size=(20,4));w=rng.uniform(.1,1,20)
    s=np.einsum('n,ni,nj->ij',w,m,m)
    out=covariance_cone_diagnostic(m,s)
    assert out['numerical_match'] and out['coefficients'].min()>=0
    np.testing.assert_allclose(out['approximation'],s,rtol=1e-8,atol=1e-9)
    growth=direction_growth(out['direction'],s)
    assert growth['expected_log_lower'] < 1e-10


def test_off_diagonal_separator_when_diagonal_powers_are_identical():
    m=np.array([[1.,1.],[-1.,-1.]])
    s=np.array([[1.,-1.],[-1.,1.]])
    np.testing.assert_array_equal(m[0]**2,np.diag(s))
    out=covariance_cone_diagnostic(m,s)
    assert not out['numerical_match'] and out['trace_signal_direction']>0
    growth=direction_growth(out['direction'],s)
    assert growth['expected_log_lower'] > .1
    norm,u=matrix_gaussian_terms(growth['weights'])
    assert np.max(np.einsum('ni,ij,nj->n',m,u,m))<1e-9
    np.testing.assert_allclose(growth['expected_log_lower'],np.sum(growth['weights']*s)+norm)


def test_moment_domain_and_external_bound_are_checked():
    with pytest.raises(ValueError):matrix_gaussian_terms(np.diag([1.,0.]))
    with pytest.raises(ValueError):matrix_log_factor([1.,2.],[3.,4.],np.eye(2)*.1,-1)
    assert matrix_log_factor([1.,2.],[3.,4.],np.zeros((2,2)),0)==0


def test_general_matrix_statistic_is_not_translation_invariant():
    # Realification of a single complex frequency: a common translation rotates
    # both exposures. Only scalar multiples of I commute with every phase.
    theta=.71;o=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    t=np.diag([.1,-.2]);x=np.array([1.,2.]);z=np.array([3.,-.5])
    assert abs(x@t@z-(o@x)@t@(o@z))>.1
    np.testing.assert_allclose(x@(.1*np.eye(2))@z,(o@x)@(.1*np.eye(2))@(o@z))


def test_near_boundary_and_singular_covariance_moments():
    t=np.diag([.999,-.998,0.]);m=np.array([.02,.05,2.])
    normalizer,u=matrix_gaussian_terms(t)
    np.testing.assert_allclose(dense_mgf(t,m,np.eye(3),np.eye(3)),
        -normalizer+m@u@m,rtol=1e-12)
    actual=dense_mgf(t,m,np.diag([0.,.5,1.]),np.diag([1.,0.,.3]))
    assert actual<=-normalizer+m@u@m


def test_off_frequency_blocks_depend_on_common_translation():
    # Two distinct frequencies and all Re coordinates followed by all Im.
    q=np.array([[1.,0.],[0.,2.]])
    phase=2*np.pi*(q@np.array([.07,.11]));c=np.diag(np.cos(phase));s=np.diag(np.sin(phase))
    shift=np.block([[c,s],[-s,c]])
    t=np.zeros((4,4));t[0,1]=t[1,0]=.1;t[2,3]=t[3,2]=.1
    x=np.array([1.,2.,3.,4.]);z=np.array([-.5,.3,1.,2.])
    assert abs(x@t@z-(shift@x)@t@(shift@z))>.1
