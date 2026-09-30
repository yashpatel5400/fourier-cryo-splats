import numpy as np
import cvxpy as cp
from fourier_splats.uq_continuous import ContinuousObservationGram
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_target_projection import target_anchored_projection,enrich_target_projection
from fourier_splats.uq_pose_exchange import solve_full_weight_cut_problem,psd_root


def test_target_projection_is_continuous_gram_psd_relaxation():
    rng=np.random.default_rng(609981);k=.8*rng.normal(size=(2,3,3));ctf=rng.normal(size=(2,3))
    gram=QuadratureObservationGram(k,ctf,.6,order=20,preconditioner_rank=0);exact=ContinuousObservationGram(k,ctf,.6)
    a,norm2=gram.target([[0,0,0]],[1],.2);m=len(a);g=np.column_stack([exact.matvec(v) for v in np.eye(m)])
    for rank in [1,5,m+1]:
        projection=target_anchored_projection(gram,a,norm2,rank=rank);f,b=projection['factor'],projection['target']
        np.testing.assert_allclose(f@b,a,rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(b@b,norm2,rtol=1e-12)
        assert np.linalg.eigvalsh(g-f@f.T).min()>-1e-8
    np.testing.assert_allclose(g,f@f.T,rtol=1e-8,atol=1e-8)


def test_full_weight_projected_conic_dual_brackets_exact_unprojected_optimum():
    rng=np.random.default_rng(609982);a=rng.normal(size=(7,4));ell=a@np.array([1.,-.3,.2,.7])+.2*rng.normal(size=7)
    q,_=np.linalg.qr(np.column_stack([ell,a[:,:1]]));f=a.T@q;b=q.T@ell
    c=np.array([[.02,.03]]);z=.1;B=1.2;scale=.3;cuts=np.zeros((0,4))
    fit=solve_full_weight_cut_problem(f,b,c,cuts,z,B,scale)
    w=cp.Variable(4)
    exact=cp.Problem(cp.Minimize(z*cp.norm(w)+B*cp.norm(ell-a@w)+c.ravel()@cp.norm(cp.vstack([w[:2],w[2:]]),axis=0)))
    exact.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9)
    candidate=fit['weights'];upper=z*np.linalg.norm(candidate)+B*np.linalg.norm(ell-a@candidate)+float(c.ravel()@np.hypot(candidate[:2],candidate[2:]))
    assert fit['full_dual_lower_bound']>0
    assert fit['full_dual_lower_bound']<=exact.value+1e-6<=upper+2e-6
    np.testing.assert_allclose(fit['full_dual_lower_bound'],fit['model_objective'],rtol=1e-5,atol=1e-6)


def test_adaptive_projection_adds_omitted_continuous_representer():
    rng=np.random.default_rng(610101);k=rng.normal(size=(3,3,3));ctf=rng.normal(size=(3,3))
    gram=QuadratureObservationGram(k,ctf,.8,order=24,preconditioner_rank=0)
    exact=ContinuousObservationGram(k,ctf,.8);a,norm2=gram.target([[0,0,0]],[1],.15)
    projection=target_anchored_projection(gram,a,norm2,rank=3)
    g=np.column_stack([exact.matvec(v) for v in np.eye(len(a))]);previous=[]
    for _ in range(5):
        w=rng.normal(size=len(a));before=np.linalg.norm(projection['target']-projection['factor'].T@w)
        diagnostic=enrich_target_projection(projection,exact,w);assert diagnostic['added']
        f,b=projection['factor'],projection['target'];previous.append(w)
        np.testing.assert_allclose(f@b,a,atol=1e-11)
        assert np.linalg.eigvalsh(g-f@f.T).min()>-1e-8
        after=np.linalg.norm(b-f.T@w)
        np.testing.assert_allclose(after**2,norm2-2*a@w+w@g@w,rtol=1e-9,atol=1e-9)
        assert after>=before-1e-10
        for old in previous:np.testing.assert_allclose((g-f@f.T)@old,0.,atol=1e-8)
    assert not enrich_target_projection(projection,exact,w)['added']
