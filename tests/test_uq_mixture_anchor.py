import numpy as np
import cvxpy as cp
from fourier_splats.uq_mixture_anchor import scaled_anchor_upper


def test_scaled_anchors_bound_independent_conic_optimum():
    rng=np.random.default_rng(2821);u=rng.uniform(.1,3.,(8,6))
    w=cp.Variable(6,nonneg=True)
    problem=cp.Problem(cp.Maximize(cp.sum(cp.log(u@w))),[cp.sum(w)==1])
    optimum=problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10)
    for anchor in [rng.normal(size=8),np.log(u.max(axis=1)),np.log(u@w.value)]:
        bound=scaled_anchor_upper(np.log(u),anchor)
        assert bound>=optimum-1e-8
        np.testing.assert_allclose(scaled_anchor_upper(np.log(u),anchor+200),bound,atol=5e-12)
        old=anchor.sum()-8+np.max(np.sum(u/np.exp(anchor[:,None]),axis=0))
        assert bound<=old+1e-10
    row_max=np.log(u.max(axis=1))
    assert scaled_anchor_upper(np.log(u),row_max)<=row_max.sum()+1e-12
    shifts=rng.normal(size=8)*100
    np.testing.assert_allclose(scaled_anchor_upper(np.log(u)+shifts[:,None],row_max+shifts),
        scaled_anchor_upper(np.log(u),row_max)+shifts.sum(),atol=1e-10)
