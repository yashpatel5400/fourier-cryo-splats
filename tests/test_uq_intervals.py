import mpmath as mp
import pytest
from fourier_splats.uq_intervals import bias_aware_half_width_stable


def test_both_tail_critical_values_against_high_precision():
    with mp.workdps(100):
        for alpha in [.05, .01, .8, 1e-10]:
            for bias in [0., .1, 1., 4., 40., 41., 100., 1e6]:
                q=bias_aware_half_width_stable(1.,bias,alpha)
                sf=lambda z:mp.erfc(z/mp.sqrt(2))/2
                failure=sf(mp.mpf(q)-mp.mpf(bias))+sf(mp.mpf(q)+mp.mpf(bias))
                assert failure<=mp.mpf(alpha), (alpha,bias,q,failure)
                assert float(mp.mpf(alpha)-failure)<1e-8
    assert bias_aware_half_width_stable(0,2)==2
    for args in [(float('nan'),1),(-1,2),(1,float('inf'))]:
        with pytest.raises(ValueError):bias_aware_half_width_stable(*args)
