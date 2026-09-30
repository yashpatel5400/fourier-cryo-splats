import numpy as np
from scipy.special import logsumexp
from fourier_splats.uq_continuous_mixture import continuous_mixture_upper
from fourier_splats.uq_mixture_refinement import SelectiveCurvatureOrbit, refine_envelope_mixture
from fourier_splats.uq_mixture_fast import matrix_mixture_fit
from fourier_splats.uq_mixture_anchor import scaled_anchor_upper


def make_orbit(zero=False):
    rng=np.random.default_rng(683219)
    plane=np.pad(rng.normal(size=(6,2)),((0,0),(0,1)))
    coeff=rng.normal(size=6)*(not zero)
    return rng,SelectiveCurvatureOrbit(plane,rng.normal(size=(3,3)),.8,coeff,
        rng.uniform(-1,1,(5,6)),rng.normal(size=(5,12)))


def test_refinement_preserves_cover_and_bounds_independent_orientation_fit():
    rng,orbit=make_orbit()
    parent=continuous_mixture_upper(orbit,initial_bins=(2,1,2),max_splits=8,
        refit_every=8,mixture_iterations=20,tolerance=1e-10)
    result=refine_envelope_mixture(orbit,parent,max_splits=48,batch_size=8,
        mixture_iterations=50,lower_iterations=100,tolerance=1e-10)
    assert result['status']=='split_limit'
    assert len(result['split_history'])==48
    low,high=result['leaf_lower'],result['leaf_upper']
    np.testing.assert_allclose(np.prod(high-low,axis=1).sum(),4*np.pi**3,rtol=1e-14)
    samples=rng.uniform(size=(1200,3))*[2*np.pi,np.pi,2*np.pi]
    kernels=[]
    for angle in samples:
        membership=np.all((angle>=low)&(angle<=high),axis=1)
        assert membership.sum()==1
        value=orbit.log_kernel(angle);kernels.append(value)
        assert np.all(value<=result['leaf_log_envelopes'][membership][0]+1e-9)
    grid=matrix_mixture_fit(np.array(kernels).T,max_iterations=100)
    assert grid['primal']+orbit.normalization<=result['log_likelihood_upper']+1e-9
    np.testing.assert_allclose(scaled_anchor_upper(result['best_cover_log_envelopes'].T,
        result['best_anchor_log_z'])+orbit.normalization,result['log_likelihood_upper'],atol=1e-9)
    finite=logsumexp(result['support_log_kernels'][result['best_feasible_ids']].T+
        np.log(result['best_feasible_weights'])[None,:],axis=1).sum()+orbit.normalization
    np.testing.assert_allclose(finite,result['log_likelihood_lower'],atol=1e-10)
    assert finite>=parent['log_likelihood_lower']-1e-9
    # Every child remains within its original parent's kernel envelope.
    for l,h,u in zip(low,high,result['leaf_log_envelopes']):
        containing=np.all((l>=parent['leaf_lower']-1e-12)&(h<=parent['leaf_upper']+1e-12),axis=1)
        assert containing.sum()==1
        assert np.all(u<=parent['leaf_log_envelopes'][containing][0]+1e-12)


def test_zero_orbit_exactness_and_unfinished_budget():
    _,orbit=make_orbit(True)
    parent=continuous_mixture_upper(orbit,initial_bins=(2,1,2),max_splits=0)
    result=refine_envelope_mixture(orbit,parent,max_splits=0)
    exact=-.5*np.sum(orbit.observed**2)+orbit.normalization
    assert result['converged']
    np.testing.assert_allclose([result['log_likelihood_lower'],result['log_likelihood_upper']],exact,atol=1e-10)
    _,orbit=make_orbit()
    parent=continuous_mixture_upper(orbit,initial_bins=(2,1,2),max_splits=0)
    result=refine_envelope_mixture(orbit,parent,max_splits=0,mixture_iterations=0,lower_iterations=0)
    assert result['log_likelihood_lower']<=result['log_likelihood_upper']
    assert not result['converged']
