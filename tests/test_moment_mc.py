import numpy as np
from scipy.stats import binom
from fourier_splats.uq_bispectrum import moment_features
from fourier_splats.uq_moment_mc import (noisy_amplitude_polynomial,cubic_interval_maximum,
    binomial_upper,rejection_count,projected_rejection_probability)


def test_noisy_amplitude_polynomial_and_supremum():
    rng=np.random.default_rng(761)
    m=rng.normal(size=(128,5))+1j*rng.normal(size=(128,5))
    e=rng.normal(size=m.shape)+1j*rng.normal(size=m.shape)
    triads=np.array([[0,1,2],[0,2,3],[1,3,4]]);w=rng.normal(size=11)
    c=noisy_amplitude_polynomial(m,e,triads,w)
    for amplitude in [-1.,.3,.9,1.,1.1,2.]:
        actual=moment_features(amplitude*m+e,triads,1.)@w
        predicted=sum(c[:,j]*amplitude**j for j in range(4))
        np.testing.assert_allclose(actual,predicted,atol=8e-14,rtol=2e-12)
    best,arg=cubic_interval_maximum(c,.9,1.1)
    dense=np.max(sum(c[:,j,None]*np.linspace(.9,1.1,10001)**j for j in range(4)),axis=1)
    np.testing.assert_allclose(best,dense,atol=2e-8,rtol=2e-9)
    np.testing.assert_allclose(best,moment_features(arg[:,None]*m+e,triads,1.)@w,atol=8e-14)


def test_cubic_degeneracy_and_cancellation():
    c=np.array([[0,0,0,0],[1,2,0,0],[0,2,-1,0],[0,0,0,1],
        [0,-3,0,1],[0,1,1e10,-1e-10]])
    best,arg=cubic_interval_maximum(c,-2.,2.)
    np.testing.assert_allclose(best[:5],[0,5,1,8,2])
    assert np.isfinite(best).all() and np.all(abs(arg)<=2)


def test_binomial_simulation_upper_and_unconditional_size():
    # Exact summation over the calibration count, including the failure event.
    m,n,alpha,delta=19,23,.05,.01
    for p in [0.,.01,.2,.5,.95,1.]:
        counts=np.arange(m+1);upper=np.array([binomial_upper(int(k),m,delta) for k in counts])
        mass=binom.pmf(counts,m,p)
        assert mass[upper<p-1e-14].sum()<=delta+1e-12
        conditional=np.array([binom.sf(rejection_count(n,u,alpha-delta)-1,n,p) for u in upper])
        assert mass@conditional<=alpha+1e-12


def test_density_ratio_domination_with_unknown_amplitude():
    # Enumerate a small orientation/noise space. Envelope allows even the
    # stronger, nonphysical noise-dependent amplitude choice.
    q=np.array([.2,.3,.5]);ratio=np.array([1.5,1.5,.5]);p=q*ratio
    assert abs(p.sum()-1)<1e-15
    events=np.array([[[0,0,1,0],[0,1,0,0]],[[0,1,0,1],[1,0,0,0]],[[0,0,0,0],[0,0,1,0]]])
    envelope=events.max(axis=1).mean(axis=1)
    assert p@envelope<=1.5*(q@envelope)
    for amplitude in [0,1]:assert p@events[:,amplitude].mean(axis=1)<=1.5*(q@envelope)


def test_binomial_critical_values_are_minimal_and_projection_is_monotone():
    for p in [0.,.003,.5,.999,1.]:
        k=rejection_count(100,p,.049)
        assert binom.sf(k-1,100,p)<=.049
        assert k==0 or binom.sf(k-2,100,p)>.049
        r=projected_rejection_probability(31,100,100,k)
        assert r['rejection_interval'][0]<=r['rejection_probability']<=r['rejection_interval'][1]
