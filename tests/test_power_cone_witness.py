import numpy as np
from fourier_splats.uq_power_cone import FinitePowerCone
from fourier_splats.uq_power_cone_witness import cone_witness,power_witness_log_upper


def test_exact_mixture_matches_powers_and_has_zero_growth_upper():
    rng=np.random.default_rng(713);p=rng.uniform(.2,2,(9,6));weights=rng.uniform(0,1,9)
    s=p.T@weights;out=cone_witness(p,s)
    assert out['equality_success']
    assert np.min(out['coefficients'])>=0
    assert out['maximum_scaled_residual']<1e-8
    assert power_witness_log_upper(s,out['approximation'],np.arange(1,7))<1e-10


def test_infeasible_target_retains_approximation_and_valid_dual():
    p=np.array([[1.,2.],[2.,4.]]);s=np.array([3.,1.]);out=cone_witness(p,s)
    assert not out['equality_success']
    assert out['maximum_scaled_residual']>.1
    upper=power_witness_log_upper(s,out['approximation'],np.ones(2))
    fitted=FinitePowerCone(2,2).solve(p,s)
    assert upper>=fitted['expected_log_lower']-1e-8
