import numpy as np
from numpy.testing import assert_allclose
from fourier_splats.uq_power_cone import FinitePowerCone, scalar_dual_maximum


def test_scalar_dual_matches_closed_form_zero_penalty():
    for s in [0., .01, 1., 10., 100.]:
        value,x=scalar_dual_maximum(s,0.)
        expected=s/(1+np.sqrt(1+s*s))
        assert_allclose(x,expected,atol=1e-13)
        assert_allclose(value,s*expected+np.log1p(-expected**2),atol=1e-12)


def test_true_mean_cone_and_rescaling_controls():
    rng=np.random.default_rng(413);powers=rng.uniform(.1,2.,(6,8));s=powers.mean(axis=0)
    fit=FinitePowerCone(6,8)
    for scaled in [powers,4*powers]:
        out=fit.solve(scaled,s)
        assert out['expected_log_lower']<1e-7
        assert out['expected_log_dual_upper']<1e-6


def test_separated_cone_and_zero_null_have_positive_growth():
    fit=FinitePowerCone(2,2);p=np.array([[1.,2.],[2.,4.]]);s=np.array([3.,1.])
    out=fit.solve(p,s)
    assert out['expected_log_lower']>.1
    assert out['gap']<1e-5
    assert np.max(p@(out['weights']/(1-out['weights'])))<1e-8
    zero=fit.solve(np.zeros((2,2)),s)
    expected=sum(scalar_dual_maximum(a,0.)[0] for a in s)
    assert_allclose(zero['expected_log_lower'],expected,atol=1e-5)
