import numpy as np
from scipy.optimize import minimize
from fourier_splats.uq_folded_ridge import folded_ridge_search,folded_width_derivatives
from fourier_splats.uq_intervals import bias_aware_half_width_stable


class DenseGram:
    def __init__(self,a,ell):
        self.a=a;self.ell=ell;self.g=a@a.T;self.shape=self.g.shape;self.diagonal=np.diag(self.g)
    def target(self,*args):return self.a@self.ell,float(self.ell@self.ell)
    def matvec(self,w):return self.g@w


def test_width_derivatives_against_finite_differences():
    for s,b in [(1.,0.),(1.,.2),(2.,3.),(.1,1.)]:
        qs,qb=folded_width_derivatives(s,b)
        step=1e-5
        ds=(bias_aware_half_width_stable(s+step,b)-bias_aware_half_width_stable(s-step,b))/(2*step)
        np.testing.assert_allclose(qs,ds,rtol=1e-6,atol=1e-7)
        if b>0:
            db=(bias_aware_half_width_stable(s,b+step)-bias_aware_half_width_stable(s,b-step))/(2*step)
            np.testing.assert_allclose(qb,db,rtol=1e-6,atol=1e-7)
        assert qs>0 and qb>=0


def test_ridge_search_brackets_independent_weight_optimization():
    rng=np.random.default_rng(381);a=rng.normal(size=(4,6));ell=rng.normal(size=6)
    gram=DenseGram(a,ell);B=2.
    result=folded_ridge_search(gram,None,None,None,B,relative_tolerance=.002,maximum_evaluations=300)
    def objective(w):return bias_aware_half_width_stable(np.linalg.norm(w),B*np.linalg.norm(ell-a.T@w))
    fit=minimize(objective,np.zeros(4),method='BFGS',options={'gtol':1e-9})
    assert result['converged']
    assert result['lower_bound']<=fit.fun+1e-8<=result['half_width']+1e-8
    assert result['half_width']/fit.fun<1.003
    q_s,q_b=folded_width_derivatives(result['noise_sd'],result['bias'])
    implied=q_s*result['bias']/(q_b*B*B*result['noise_sd'])
    assert abs(np.log(implied/result['ridge']))<.1


def test_bounds_cover_a_dense_path_including_extreme_ridges():
    rng=np.random.default_rng(612);a=rng.normal(size=(3,4));ell=rng.normal(size=4)
    result=folded_ridge_search(DenseGram(a,ell),None,None,None,1.3,maximum_evaluations=20)
    values=[]
    for ridge in np.geomspace(1e-8,1e8,301):
        w=np.linalg.solve(a@a.T+ridge*np.eye(3),a@ell)
        values.append(bias_aware_half_width_stable(np.linalg.norm(w),1.3*np.linalg.norm(ell-a.T@w)))
    assert result['global_lower_bound']<=min(values)+1e-8
    assert result['half_width']<=result['no_data_half_width']


def test_unobserved_target_selects_no_data_endpoint():
    gram=DenseGram(np.array([[1.,0.]]),np.array([0.,1.]))
    result=folded_ridge_search(gram,None,None,None,2.,maximum_evaluations=12)
    assert result['half_width']==2. and result['converged']
    np.testing.assert_array_equal(result['weights'],[0.])
