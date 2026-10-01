"""Exact finite-space verification of the retrospective reciprocal premise."""
import itertools
import numpy as np


def expectation(f,q,planted,n):
    result=0.
    for indices in itertools.product(range(len(q)),repeat=n):
        probability=planted[indices[0]]*np.prod([q[k] for k in indices[1:]])
        estimate=sum(f[k]/q[k] for k in indices)/n
        result+=probability/estimate
    return result


def test_posterior_replacement_is_reciprocal_unbiased_exactly():
    f=np.array([.013,1.2,3.7]);q=np.array([.85,.1,.05])
    for n in [1,2,4]:
        np.testing.assert_allclose(expectation(f,q,f/f.sum(),n),1/f.sum(),rtol=2e-15,atol=0)


def test_wrong_state_and_fixed_truth_do_not_supply_the_identity():
    f=np.array([.013,1.2,3.7]);q=np.array([.85,.1,.05])
    wrong_state=np.array([.8,.15,.05])
    fixed_pose=np.array([0.,0.,1.])
    assert abs(expectation(f,q,wrong_state,4)-1/f.sum())>.1
    assert abs(expectation(f,q,fixed_pose,4)-1/f.sum())>.1
