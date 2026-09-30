import numpy as np
import cvxpy as cp
from scipy.special import logsumexp
from fourier_splats.uq_mixture_fast import matrix_mixture_fit
from fourier_splats.uq_mixture_anchor import scaled_anchor_upper


def test_matrix_em_certificate_and_warm_start_against_conic_fit():
    rng=np.random.default_rng(63717)
    u=rng.uniform(.05,2,(25,12));w=cp.Variable(12,nonneg=True)
    problem=cp.Problem(cp.Maximize(cp.sum(cp.log(u@w))),[cp.sum(w)==1])
    exact=problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10)
    first=matrix_mixture_fit(np.log(u),max_iterations=10,tolerance=1e-8)
    final=matrix_mixture_fit(np.log(u),weights=first['weights'],max_iterations=5000,tolerance=1e-5)
    for result in [first,final]:
        assert result['primal']<=exact+1e-7<=result['upper']+2e-7
        np.testing.assert_allclose(logsumexp(np.log(u)+np.log(result['weights'])[None,:],axis=1).sum(),result['primal'],atol=1e-12)
        np.testing.assert_allclose(scaled_anchor_upper(np.log(u),result['upper_log_anchor']),result['upper'],atol=1e-10)
    assert final['converged']
    assert final['primal']>=first['primal']-1e-12


def test_original_extreme_log_kernels_and_zero_iteration_bounds():
    values=np.array([[0.,-3000.,-1400.],[-1000.,0.,-2000.],[-1400.,-5000.,0.]])
    fit=matrix_mixture_fit(values,max_iterations=0)
    exact=-3*np.log(3)
    np.testing.assert_allclose(fit['primal'],exact,atol=1e-12)
    np.testing.assert_allclose(fit['upper'],exact,atol=1e-12)
    shifts=np.array([-100000.,2000.,100.])
    moved=matrix_mixture_fit(values+shifts[:,None],max_iterations=0)
    np.testing.assert_allclose(moved['upper'],fit['upper']+shifts.sum(),atol=1e-9)
