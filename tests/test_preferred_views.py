import numpy as np
import pytest
from fourier_splats.uq_preferred_views import sample_cap_views


@pytest.mark.parametrize('axis,kappa',[(0,1.1),(1,2.),(2,5.)])
def test_known_viewing_density_and_moments(axis,kappa):
    r,info=sample_cap_views(20000,axis,kappa,np.random.default_rng(100+axis))
    cut=1-1/kappa;coordinate=r[:,2,axis]
    assert np.min(abs(coordinate))>=cut
    np.testing.assert_allclose(r@r.transpose(0,2,1),np.broadcast_to(np.eye(3),r.shape),atol=2e-15)
    np.testing.assert_allclose(np.linalg.det(r),1,atol=2e-15)
    assert abs(coordinate.mean())<.025
    assert abs(np.mean(coordinate**2)-(1+cut+cut*cut)/3)<.01
    assert abs(info['accepted']/info['proposed']-1/kappa)<.012
    # In-plane rotations remain uniform: cap restriction involves only row 2.
    assert np.max(abs(r[:,:2,:].mean(axis=0)))<.025
    replay,_=sample_cap_views(20000,axis,kappa,np.random.default_rng(100+axis))
    np.testing.assert_array_equal(r,replay)
