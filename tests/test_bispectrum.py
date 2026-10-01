import itertools
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize
from fourier_splats.uq_bispectrum import (frequency_triads, moment_features,
    MomentContrast, project_simplex, moment_hull_projection, amplitude_maximum)


def test_continuous_translation_invariance():
    q=np.array([[1,0],[0,1],[1,1],[2,1],[1,2]],float)
    triads=frequency_triads(q,maximum=None)
    assert len(triads)==3
    rng=np.random.default_rng(971)
    values=rng.normal(size=(11,len(q)))+1j*rng.normal(size=(11,len(q)))
    shifts=rng.uniform(-100,100,(11,2))
    shifted=values*np.exp(-2j*np.pi*(shifts@q.T))
    np.testing.assert_allclose(moment_features(shifted,triads),moment_features(values,triads),atol=1e-12)


def test_gaussian_variance_against_exact_six_dimensional_quadrature():
    triads=np.array([[0,1,2]])
    weights=np.array([.3,-.1,.9,.4,-.7])
    means=np.array([[.5+.8j,-.2+.3j,1.2-.9j]])
    contrast=MomentContrast(3,triads,weights)
    x,w=hermgauss(4)
    indices=np.array(list(itertools.product(range(4),repeat=6)))
    nodes=np.sqrt(2)*x[indices];mass=np.prod(w[indices]/np.sqrt(np.pi),axis=1)
    for v in [0.,.2,1.,2.]:
        noise=np.sqrt(v)*(nodes[:,:3]+1j*nodes[:,3:])
        scores=moment_features(means+noise,triads,noise_variance=v)@weights
        expected,variance=contrast.mean_variance(means,v)
        np.testing.assert_allclose(mass@scores,expected[0],atol=2e-13)
        np.testing.assert_allclose(mass@((scores-expected[0])**2),variance[0],rtol=2e-13,atol=2e-13)


def test_shared_frequency_covariances_against_monte_carlo():
    q=np.array([[1,0],[0,1],[1,1],[2,1],[1,2]],float)
    t=frequency_triads(q,maximum=None)
    rng=np.random.default_rng(9712);mean=rng.normal(size=5)+1j*rng.normal(size=5)
    w=rng.normal(size=5+2*len(t));c=MomentContrast(5,t,w)
    noise=rng.normal(size=(200000,5))+1j*rng.normal(size=(200000,5))
    scores=moment_features(mean+noise,t,noise_variance=1.)@w
    expected,variance=c.mean_variance(mean)
    assert abs(scores.mean()-expected[0]) < 5*np.sqrt(variance[0]/len(scores))
    np.testing.assert_allclose(scores.var(),variance[0],rtol=.025)


def test_hull_bracket_against_independent_quadratic_solver():
    rng=np.random.default_rng(771)
    atoms=rng.normal(size=(14,5));target=np.array([3.,-2,1,0,.5])
    result=moment_hull_projection(atoms,target,max_iterations=6000,maximum_seconds=10.,distance_tolerance=1e-8)
    objective=lambda w:.5*np.sum((w@atoms-target)**2)
    fit=minimize(objective,np.full(14,1/14),jac=lambda w:atoms@(w@atoms-target),
        method='SLSQP',bounds=[(0,1)]*14,constraints=[{'type':'eq','fun':lambda w:w.sum()-1,
        'jac':lambda w:np.ones(14)}],options={'ftol':1e-13,'maxiter':1000})
    assert fit.success
    optimum=np.sqrt(2*fit.fun)
    assert result['distance_lower'] <= optimum+1e-8 <= result['distance_upper']+2e-8
    assert result['gap']<1e-7
    assert result['weights'].min()>=0
    np.testing.assert_allclose(result['weights'].sum(),1.,atol=1e-14)


def test_hull_membership_and_simplex():
    atoms=np.array([[0.,0.],[1.,0.],[0.,1.]])
    r=moment_hull_projection(atoms,np.array([.2,.3]),distance_tolerance=1e-8)
    assert r['distance_upper']<1e-8 and r['distance_lower']<1e-8
    np.testing.assert_allclose(project_simplex(np.array([3.,-4.,2.])),[1.,0.,0.])


def test_amplitude_polynomial_and_global_profile():
    rng=np.random.default_rng(551)
    means=rng.normal(size=(7,3))+1j*rng.normal(size=(7,3))
    c=MomentContrast(3,[[0,1,2]],rng.normal(size=5))
    polynomial=c.variance_polynomial(means)
    for amplitude in [.2,.9,1.,1.1,2.]:
        actual=c.mean_variance(amplitude*means)[1]
        expected=polynomial@np.array([amplitude**j for j in range(5)])
        np.testing.assert_allclose(actual,expected,rtol=1e-13)
    examples=np.array([[0,0,2,-1,0],[1,0,-3,2,-.3],[0,0,0,0,0]],float)
    maximum,arg=amplitude_maximum(examples,0,4)
    grid=np.linspace(0,4,10001)
    values=examples@np.stack([grid**j for j in range(5)])
    assert np.all(maximum>=values.max(axis=1)-1e-12)
    assert np.max(maximum-values.max(axis=1))<1e-6
    np.testing.assert_allclose(arg[0],4/3,atol=1e-14)
