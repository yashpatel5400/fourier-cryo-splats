import numpy as np
from fourier_splats.uq_shared_density import grouped_spectral_bound,shared_density_bounds


def test_spectral_bound_retains_orthogonality_and_covers_coupled_quadratics():
    # Orthogonal block images make the sum-of-norms bound unnecessarily large.
    eye=np.eye(4);blocks=[eye[:,i:i+1] for i in range(4)]
    np.testing.assert_allclose(grouped_spectral_bound(blocks),2,rtol=1e-10)
    assert grouped_spectral_bound([np.zeros((4,2))])==0
    rng=np.random.default_rng(609355);p,n,d=11,4,3
    h=rng.normal(size=p);linear=[rng.normal(size=(p,d)) for _ in range(n)]
    quadratic=[rng.normal(size=(p,d*d)) for _ in range(n)]
    bounds=shared_density_bounds(h,linear,quadratic)
    assert bounds['minimum_valid_bound']<=bounds['block_triangle']
    for _ in range(100):
        u=rng.normal(size=(n,d));u/=np.maximum(1,np.linalg.norm(u,axis=1,keepdims=True))
        value=h+sum(a@v+b@np.outer(v,v).ravel() for a,b,v in zip(linear,quadratic,u))
        assert np.linalg.norm(value)<=min(bounds.values())+1e-10
