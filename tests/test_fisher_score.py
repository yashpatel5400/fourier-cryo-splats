import numpy as np
from scipy.linalg import null_space
from fourier_splats.uq_fisher_score import constrained_fisher_direction


def test_constrained_discriminant_matches_independent_nullspace_solution():
    rng=np.random.default_rng(82);p=13;u=rng.normal(size=p);v=u+rng.normal(size=p)
    a=rng.normal(size=(50,p));cov=a.T@a/50
    got,info=constrained_fisher_direction(u,v,cov,5)
    constraints=np.zeros((p,2));constraints[:5,0]=u[:5];constraints[5:,1]=u[5:]
    basis=null_space(constraints.T);reg=.9*cov+.1*np.trace(cov)/p*np.eye(p)
    expected=basis@np.linalg.solve(basis.T@reg@basis,basis.T@(v-u));expected/=np.linalg.norm(expected)
    np.testing.assert_allclose(got,expected,atol=1e-12)
    assert np.max(abs(constraints.T@got))<1e-12
    for _ in range(100):
        z=basis@rng.normal(size=p-2)
        assert (z@(v-u))**2/(z@reg@z)<=info['regularized_squared_separation']+1e-12


def test_identity_shrinkage_reduces_to_euclidean_nuisance_projection():
    rng=np.random.default_rng(2);u=rng.normal(size=8);d=rng.normal(size=8)
    got,_=constrained_fisher_direction(u,u+d,np.diag(np.arange(1,9)),3,shrinkage=1.)
    for sl in [slice(0,3),slice(3,8)]:d[sl]-=np.dot(d[sl],u[sl])/np.dot(u[sl],u[sl])*u[sl]
    np.testing.assert_allclose(got,d/np.linalg.norm(d),atol=1e-12)


def test_amplitude_only_alternative_is_not_identifiable_after_projection():
    u=np.array([1.,2.,3.,4.,5.]);v=np.r_[.9**2*u[:3],.9**3*u[3:]]
    h,r=constrained_fisher_direction(u,v,np.eye(5),3)
    assert r['zero_direction'] and np.all(h==0)
