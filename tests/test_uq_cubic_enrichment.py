import numpy as np
from test_uq_cubic_design_review import two_particle_fixture, conic_with_explicit_particle_packing
from fourier_splats.uq_cubic_subspace import orthonormal_basis
from fourier_splats.uq_cubic_enrichment import enrich_basis, optimize_enriched_cubic


def test_enrichment_keeps_span_and_rejects_dependent_directions():
    rng=np.random.default_rng(952101); q,_=orthonormal_basis(rng.normal(size=(17,3)))
    extended,records=enrich_basis(q,[q[:,0],rng.normal(size=17),np.zeros(17)])
    assert [r['added'] for r in records]==[False,True,False]
    np.testing.assert_allclose(extended[:,:3],q,atol=0)
    np.testing.assert_allclose(extended.T@extended,np.eye(4),atol=1e-12)


def test_enrichment_against_independent_full_dense_conic_problem():
    obj,w=two_particle_fixture(); obj.L*=.02
    obj.remainder.multipliers*=.02; obj.remainder.residual*=.02
    reference=conic_with_explicit_particle_packing(obj, solver='SCS')
    # The independent conic helper returns its objective and optimizer; verify
    # this richer fixture actually has a nonzero optimum.
    expected, solution, status = reference
    assert status == 'optimal' and np.linalg.norm(solution) > 1e-4
    q,_=orthonormal_basis(w[:,None])
    result=optimize_enriched_cubic(obj,q,w,optimization_seed=952102,certificate_seed=952103,
        outer_rounds=8,evaluations_per_round=15,design_seconds=300,separation_tolerance=1e-7,
        modes=2,dense_oracle=True,power_iterations=50)
    assert result['selected_approximate_objective'] < result['optimization_history'][0]['guide_objective']
    assert abs(result['selected_approximate_objective']-expected) < 5e-4
    assert result['enriched_basis'].shape[1]>1
    assert result['sum_objective_upper']>=result['dual_lower_bound']
    assert result['spectral_upper_certificate']['seed']==952103
    print('ENRICHMENT_CHECK',expected,result['selected_approximate_objective'],result['enriched_basis'].shape)
