import itertools
import numpy as np
from fourier_splats.uq_moment_mc import cubic_interval_maximum
from fourier_splats.uq_view_variance import (cubic_cell_events,grouped_amplitude_envelope,
    empirical_bernstein_radius,viewing_probability_bounds)


def test_cell_events_against_independent_polynomial_roots():
    rng=np.random.default_rng(934);c=rng.normal(size=(128,4))
    c[:4]=[[1,0,0,0],[0,1,0,0],[0,2,-1,0],[0,0,0,1]]
    events=cubic_cell_events(c,-1.,2.,16,.7);edges=np.linspace(-1,2,17)
    for i,coefficients in enumerate(c):
        p=np.polynomial.Polynomial(coefficients);roots=p.deriv().roots()
        for j,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
            points=[a,b]+[r.real for r in roots if abs(r.imag)<1e-10 and a<=r.real<=b]
            assert events[i,j]==(np.max(p(points))>.7)


def test_grouping_order_and_continuous_amplitude_domination():
    rng=np.random.default_rng(2);c=rng.normal(size=(7,2,16,4))
    grouped,individual=grouped_amplitude_envelope(c,.9,1.1,8,.3)
    assert np.all(grouped<=individual)
    exact=cubic_interval_maximum(c.reshape(-1,4),.9,1.1)[0].reshape(7,2,16)
    np.testing.assert_array_equal(individual,(exact>.3).mean(axis=2))
    for a in np.linspace(.9,1.1,21):
        score=sum(c[...,j]*a**j for j in range(4))
        assert np.all((score>.3).mean(axis=2)<=grouped)


def test_nested_envelope_expectation_and_conditional_product_identity():
    # Enumerate two noise draws per group on a four-state noise space.
    probabilities=np.array([.1,.2,.3,.4]);events=np.array([[1,0],[0,1],[0,0],[1,1]])
    group=[];mass=[]
    for a,b in itertools.product(range(4),repeat=2):
        group.append(events[[a,b]].mean(axis=0).max());mass.append(probabilities[a]*probabilities[b])
    group=np.array(group);mass=np.array(mass);eta=group@mass
    assert eta>=np.max(probabilities@events)
    center=.4
    product=np.outer(group-center,group-center)
    np.testing.assert_allclose(mass@product@mass,(eta-center)**2)


def test_density_ratio_variance_bound_exact_finite_orientation_space():
    rng=np.random.default_rng(55)
    for _ in range(50):
        q=rng.dirichlet(np.ones(8));eta=rng.uniform(0,1,8);w=rng.uniform(.2,2,8);w/=q@w
        kappa=w.max();mu=q@eta;variance=q@(eta-mu)**2
        assert q@(w-1)**2<=kappa-1+1e-12
        assert q@(w*eta)<=mu+np.sqrt((kappa-1)*variance)+1e-12


def test_bernstein_scaling_and_variance_bound_formula():
    rng=np.random.default_rng(8);x=rng.uniform(0,1,100)
    radius=empirical_bernstein_radius(x,0,1,.001)
    np.testing.assert_allclose(empirical_bernstein_radius(3*x-2,-2,1,.001),3*radius)
    # Constant conditional probability represented by independent noise groups.
    groups=rng.binomial(32,.4,size=(20000,2))/32
    bound=viewing_probability_bounds(groups,groups,.4,[1,1.1,2])
    assert bound['view_variance_upper']>=0
    assert bound['joint_mean_interval'][0]<.4<bound['joint_mean_interval'][1]
    for row in bound['bounds']:
        assert .4<=row['view_variance']<=1
    assert bound['bounds'][1]['view_variance']<bound['bounds'][1]['grouped_ratio']
