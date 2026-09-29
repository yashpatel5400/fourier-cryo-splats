import numpy as np
import pytest
from scipy.stats import norm
from fourier_splats.uq_subspace import enrich_certificate,matrix_free_certificate


class DenseOperator:
    def __init__(self,a):self.a=a;self.shape=a.shape
    def forward(self,v):return self.a@v
    def adjoint(self,w):return self.a.T@w
    def forward_columns(self,s):return self.a@s


def test_enrichment_against_ambient_conic_optimum_and_adversary():
    cp=pytest.importorskip('cvxpy');rng=np.random.default_rng(909)
    a=rng.normal(size=(35,25));ell=rng.normal(size=25);B=2.
    initial=rng.normal(size=(25,4));op=DenseOperator(a)
    fit=enrich_certificate(op,ell,initial,B,max_rounds=30,rtol=1e-4,inner_rtol=1e-5,inner_maxiter=400)
    w=cp.Variable(35);cost=norm.ppf(.975)*cp.norm(w)+B*cp.norm(ell-a.T@w)
    optimum=cp.Problem(cp.Minimize(cost)).solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10)
    assert fit.converged
    assert fit.dual_lower_bound<=optimum+1e-7
    assert fit.objective>=optimum-1e-7
    assert fit.objective-optimum<1e-4*optimum
    assert np.all(np.diff([r['best_ambient_objective'] for r in fit.history])<=1e-12)
    residual=ell-a.T@fit.weights;delta=B*residual/np.linalg.norm(residual)
    bias=fit.weights@a@delta-ell@delta
    actual=norm.cdf((fit.half_width-bias)/fit.noise_sd)-norm.cdf((-fit.half_width-bias)/fit.noise_sd)
    assert actual==pytest.approx(.95,abs=1e-8)
    assert fit.history[0]['ambient_half_width']>fit.history[0]['restricted_half_width']
    full=matrix_free_certificate(op,ell,B,maxiter=100,rtol=1e-4)
    assert full.converged
    assert full.dual_lower_bound<=optimum+1e-7
    assert full.objective-optimum<1e-4*optimum
