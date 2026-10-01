import itertools
import numpy as np
from scipy.optimize import linprog
from fourier_splats.uq_view_risk import empirical_tail_bound,risk_baselines,paired_variance_decomposition
from fourier_splats.uq_view_variance import grouped_amplitude_envelope,cubic_cell_events,viewing_probability_bounds


def test_cvar_against_fractional_density_linear_program():
    x=np.array([.1,.1,.2,.5,.9,.9,1.])
    for kappa in [1.,1.01,1.1,2.,5.,100.]:
        fit=linprog(-x/len(x),A_eq=np.ones((1,len(x)))/len(x),b_eq=[1.],bounds=[(0,kappa)]*len(x),method='highs')
        assert fit.success
        np.testing.assert_allclose(empirical_tail_bound(x,kappa)[0],-fit.fun,atol=1e-12)


def test_dkw_objective_independent_breakpoint_enumeration():
    rng=np.random.default_rng(10);x=rng.uniform(0,1,80)
    for kappa in [1.,1.1,2.,5.]:
        epsilon=.019
        exact=min(t+kappa*(np.maximum(x-t,0).mean()+epsilon*(1-t)) for t in np.r_[0.,x,1.])
        np.testing.assert_allclose(empirical_tail_bound(x,kappa,epsilon)[0],exact,atol=2e-15)


def test_conditional_jensen_for_tail_functional():
    conditional=np.array([[0.,.1,.2],[.1,.3,.5],[.4,.7,1.]])
    eta=conditional.mean(axis=1)
    for kappa in [1.,1.1,2.,5.]:
        assert empirical_tail_bound(eta,kappa)[0]<=empirical_tail_bound(conditional.ravel(),kappa)[0]+1e-15


def test_noise_adaptive_amplitude_is_not_grouped_null():
    # f(a+e)=(a+e-1)^2, e uniform on {-1,+1}, threshold 2.25.
    # A chosen from e makes the event certain. Any fixed amplitude gives <=1/2.
    noises=np.array(list(itertools.product([-1.,1.],repeat=2)))
    c=np.zeros((len(noises),2,2,4));c[...,0]=(noises[:,None,:]-1)**2
    c[...,1]=2*(noises[:,None,:]-1);c[...,2]=1
    grouped,individual=grouped_amplitude_envelope(c,0,2,4,2.25)
    assert np.mean(grouped)==.75 and np.mean(individual)==1.
    assert np.max([np.mean((a+np.array([-1.,1.])-1)**2>2.25) for a in np.linspace(0,2,101)])==.5


def test_adversarial_cubic_events():
    # Quadratic peak at the exact shared edge; tiny cubic; double stationary root.
    coefficients=np.array([[0.,2.,-1.,0.],[0.,2.,-1.,1e-300],[0.,3.,-3.,1.]])
    events=cubic_cell_events(coefficients,0,2,4,.999999)
    assert events[0].tolist()==[False,True,True,False]
    np.testing.assert_array_equal(events[0],events[1])
    assert events[2].tolist()==[False,True,True,True]


def test_variance_decomposition_replays_production():
    rng=np.random.default_rng(20);x=rng.binomial(32,.4,size=(3000,2))/32
    production=viewing_probability_bounds(x,x,.4,[1,1.1,2])
    decomposition=paired_variance_decomposition(x,.4)
    np.testing.assert_allclose(decomposition['variance_upper'],production['view_variance_upper'],atol=1e-16)
    baseline=risk_baselines(x,.4,[1,1.1,2])
    for row in baseline['bounds']:
        for key in ['cvar_dkw','cvar_split','mean_only','unpaired_variance']:
            assert .4<=row[key]<=1.
