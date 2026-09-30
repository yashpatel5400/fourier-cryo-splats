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


def test_full_weight_nonempty_spectral_cuts_bracket_unprojected_optimum():
    rng=np.random.default_rng(610201);m=4;a=rng.normal(size=(8,m))
    ell=a@np.array([1.,-.3,.2,.7])+.1*rng.normal(size=8)
    q,_=np.linalg.qr(np.column_stack([ell,a[:,:2]]));f=a.T@q;b=q.T@ell
    assert f.shape[1]<m+1
    fields=.4*rng.normal(size=(m,5,3));c=np.array([[.03,.05]]);z=.2;B=1.2;scale=.5
    w=cp.Variable(m);matrix=sum(w[j]*fields[j] for j in range(m))
    exact=cp.Problem(cp.Minimize(z*cp.norm(w)+B*cp.norm(ell-a@w)+scale*cp.norm(matrix,2)+c.ravel()@cp.norm(cp.vstack([w[:2],w[2:]]),axis=0)))
    exact.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9)
    def support(candidate):
        u,_,vh=np.linalg.svd(np.einsum('j,jab->ab',candidate,fields),full_matrices=False)
        return np.einsum('a,jab,b->j',u[:,0],fields,vh[0])
    initial=support(np.ones(m));cuts=[initial,-initial]
    for _ in range(15):
        fit=solve_full_weight_cut_problem(f,b,c,np.asarray(cuts),z,B,scale)
        candidate=fit['weights'];field=np.einsum('j,jab->ab',candidate,fields)
        upper=z*np.linalg.norm(candidate)+B*np.linalg.norm(ell-a@candidate)+scale*np.linalg.norm(field,2)+c.ravel()@np.hypot(candidate[:2],candidate[2:])
        assert 0<fit['full_dual_lower_bound']<=exact.value+2e-6<=upper+4e-6
        assert np.linalg.norm(fit['density_dual'])<=B*(1+1e-12)
        defect=f@fit['density_dual']-fit['pose_support']-fit['cubic_support']
        assert fit['dual_scale']*np.linalg.norm(defect)<=z*(1+1e-12)
        for other in rng.normal(size=(5,m)):
            assert fit['pose_support']@other<=scale*np.linalg.norm(np.einsum('j,jab->ab',other,fields),2)+1e-10
            assert fit['cubic_support']@other<=c.ravel()@np.hypot(other[:2],other[2:])+1e-10
        new=support(candidate);cuts.extend([new,-new])
