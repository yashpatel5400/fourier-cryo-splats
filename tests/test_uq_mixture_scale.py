import numpy as np
import pytest
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from fourier_splats.uq_mixture_scale import profile_common_scale


def check_partition(result):
    intervals = sorted(result['nodes'][i]['scale_interval'] for i in result['leaves'])
    assert intervals[0][0] == result['bracket'][0]
    assert intervals[-1][1] == result['bracket'][1]
    assert all(a[1] == b[0] for a,b in zip(intervals,intervals[1:]))


def test_one_component_common_scale_against_exact_mle():
    d = 7; r = np.array([1.2,2.8,3.,2.1,4.5])[:,None]
    scale = np.sqrt(np.mean(r*r)/d)
    exact = float(np.sum(-d*np.log(scale*np.sqrt(2*np.pi))-.5*(r/scale)**2))
    fit = profile_common_scale(r,d,tolerance=.005,max_nodes=129)
    assert fit['converged']
    assert fit['feasible_log_likelihood'] <= exact+1e-12
    assert exact <= fit['log_likelihood_upper']+1e-12
    assert fit['log_likelihood_upper']-exact < .005
    check_partition(fit)


def test_two_component_bounds_contain_independent_scale_grid():
    rng = np.random.default_rng(6371);r = rng.uniform(.5,3.,size=(9,2));d = 3
    fit = profile_common_scale(r,d,tolerance=.02,max_nodes=129,
        mixture_iterations=5000,mixture_tolerance=1e-7)
    points=[]
    # Include scales far outside the data-dependent bracket: tails must not
    # exceed its global upper bound either.
    for scale in np.geomspace(fit['bracket'][0]/5,fit['bracket'][1]*5,401):
        log_kernel = -d*np.log(scale*np.sqrt(2*np.pi))-.5*(r/scale)**2
        def objective(w):
            weights=np.array([w,1-w])
            with np.errstate(divide='ignore'):
                return -float(logsumexp(log_kernel+np.log(weights)[None,:],axis=1).sum())
        result=minimize_scalar(objective,bounds=(0,1),method='bounded',options={'xatol':1e-12})
        points.append(max(-result.fun,-objective(0.),-objective(1.)))
    assert max(points) <= fit['log_likelihood_upper']+1e-8
    assert fit['log_likelihood_upper']-max(points) < .03
    check_partition(fit)
    limited=profile_common_scale(r,d,max_nodes=1,tolerance=1e-12)
    assert not limited['converged'] and limited['status']=='node_limit'
    assert max(points) <= limited['log_likelihood_upper']
    degenerate=profile_common_scale(np.array([[0.,1.]]),d)
    assert degenerate['status']=='zero_distance_unresolved'
    assert np.isposinf(degenerate['log_likelihood_upper'])
