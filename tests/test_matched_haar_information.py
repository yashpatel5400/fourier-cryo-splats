"""Independent small-mixture checks for the information diagnostic."""
import importlib.util
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

spec=importlib.util.spec_from_file_location('haar_diagnostic',Path(__file__).resolve().parents[1]/'scripts/analyze_matched_haar_likelihood.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_prefix_diagnostics_and_extreme_offsets():
    rng=np.random.default_rng(934);x=rng.normal(size=(3,19))*5
    x[0,:7]-=10000;x[0,7:]+=10000
    levels=[3,7,11,19];logs,ess,weights=m.mixture_diagnostics(x,levels)
    for i,n in enumerate(levels):
        lse=logsumexp(x[:,:n],axis=1);w=np.exp(x[:,:n]-lse[:,None])
        np.testing.assert_allclose(logs[i],lse-np.log(n),atol=3e-12,rtol=0)
        np.testing.assert_allclose(ess[i],1/np.sum(w*w,axis=1),atol=3e-10,rtol=0)
        np.testing.assert_allclose(weights[i],w.max(axis=1),atol=3e-12,rtol=0)


def test_two_state_mixture_ratios_match_normal_density():
    rng=np.random.default_rng(285);y=rng.normal(size=(5,4));r=rng.normal(size=(5,4));m0=rng.normal(size=(9,4));mr=rng.normal(size=(9,4))
    for delta in [0.,.25,1.]:
        blocks=m.kernel_blocks(y,r,m0,mr,delta)
        for state in [0,1]:
            yy=y-state*delta*r
            direct=[]
            for means in [m0,m0-delta*mr]:
                matrix=np.array([multivariate_normal.logpdf(yy,mean=mu,cov=np.eye(4)) for mu in means]).T
                direct.append(logsumexp(matrix,axis=1)-np.log(len(means)))
            our=logsumexp(blocks[2*state+1],axis=1)-logsumexp(blocks[2*state],axis=1)
            np.testing.assert_allclose(our,direct[1]-direct[0],atol=3e-14,rtol=0)
