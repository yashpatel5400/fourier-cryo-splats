import numpy as np
from fourier_splats.uq_baselines import gaussian_reference,prior_predictive_coverage,unregularized_reference


def test_gaussian_posterior_has_correct_prior_predictive_coverage():
    rng=np.random.default_rng(615)
    a=rng.normal(size=(35,9));ell=rng.normal(size=9);j=rng.normal(size=(5,7,2))
    for pose in [0.,.5,3.]:
        baseline=gaussian_reference(a,ell,.4,j,pose)
        coverage=prior_predictive_coverage(a,ell,.4,baseline,j,pose)
        assert abs(coverage-.95)<1e-12


def test_unregularized_does_not_claim_identified_nullspace_functional():
    a=np.array([[1.,0],[2.,0],[3.,0]])
    valid=unregularized_reference(a,[1.,0]);invalid=unregularized_reference(a,[0.,1.])
    assert valid['identifiable'] and valid['half_width']>0
    assert not invalid['identifiable'] and invalid['half_width'] is None
