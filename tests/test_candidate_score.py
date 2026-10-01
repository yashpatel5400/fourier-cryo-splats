import numpy as np
from fourier_splats.uq_candidate_score import candidate_region,candidate_moment_direction
from fourier_splats.uq_bispectrum import moment_features


def test_candidate_peak_requires_no_reference_and_is_deterministic():
    r=np.zeros((32,32,32));r[19,16,12]=1
    mask,info=candidate_region(r,320,10)
    assert info['center_index_zyx']==[19,16,12]
    assert mask[19,16,12]==1 and np.all((mask>0)&(mask<=1))
    np.testing.assert_array_equal(candidate_region(r,320,10)[0],mask)


def test_scale_orthogonality_holds_for_all_global_amplitudes_on_training_law():
    rng=np.random.default_rng(42);m=rng.normal(size=(300,4))+1j*rng.normal(size=(300,4))
    g=.3*m+.1*(rng.normal(size=m.shape)+1j*rng.normal(size=m.shape));triads=np.array([[0,1,2]])
    direction,info=candidate_moment_direction(m,g,triads,True)
    assert not info['zero_direction'] and info['training_design_alternative_mean']>info['training_null_mean']
    for a in [.5,.9,1,1.1,2]:
        assert abs(moment_features(a*m,triads).mean(axis=0)@direction)<1e-12


def test_global_scale_change_is_removed_by_orthogonalization():
    rng=np.random.default_rng(4);m=rng.normal(size=(100,4))+1j*rng.normal(size=(100,4));triads=np.array([[0,1,2]])
    direction,info=candidate_moment_direction(m,.2*m,triads,True)
    assert info['zero_direction'] and np.all(direction==0)
