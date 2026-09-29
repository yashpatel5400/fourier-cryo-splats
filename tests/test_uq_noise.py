import numpy as np
from scipy.stats import chi2
from fourier_splats.uq_noise import gaussian_scale_upper


def test_scale_upper_tail_calibration_and_noncentral_signal():
    rng=np.random.default_rng(9743);sigma=2.3;beta=.05
    noise=rng.normal(scale=sigma,size=(100000,32))
    for mean in [0.,sigma/2]:
        upper=gaussian_scale_upper(noise+mean,beta)
        miss=np.mean(upper<sigma)
        if mean==0:assert abs(miss-beta)<.002
        else:assert miss<beta
    # An exact constructed calibration energy checks the confidence convention.
    samples=np.zeros(32);samples[0]=sigma*np.sqrt(chi2.ppf(beta,32))
    np.testing.assert_allclose(gaussian_scale_upper(samples,beta),sigma,rtol=1e-14)
