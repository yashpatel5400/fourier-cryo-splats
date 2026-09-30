import numpy as np
import cvxpy as cp
from fourier_splats.uq_pose_exchange import solve_pose_cut_problem,psd_root


def test_spectral_cut_exchange_matches_independent_full_conic_solution():
    rng=np.random.default_rng(609931);m=4;basis=np.eye(m)
    a=rng.normal(size=(8,m));ell=a@np.array([1.,-.3,.2,.7])+.05*rng.normal(size=8);joint=np.column_stack([ell,-a]);gram=joint.T@joint
    fields=rng.normal(size=(m,5,3));c=np.array([[.11,.17]]);z=.6;B=1.2;scale=.7
    w=cp.Variable(m);field=sum(w[j]*fields[j] for j in range(m))
    exact=cp.Problem(cp.Minimize(z*cp.norm(w)+B*cp.norm(ell-a@w)+scale*cp.norm(field,2)+c.ravel()@cp.norm(cp.vstack([w[:2],w[2:]]),axis=0)))
    exact.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9)
    assert np.linalg.norm(w.value)>.01
    cuts=[];history=[]
    for _ in range(35):
        fit=solve_pose_cut_problem(basis,gram,c,np.asarray(cuts).reshape(-1,m),z,B,scale)
        candidate=fit['weights'];f=np.einsum('j,jab->ab',candidate,fields)
        u,s,vh=np.linalg.svd(f,full_matrices=False)
        support=np.einsum('a,jab,b->j',u[:,0],fields,vh[0])
        objective=z*np.linalg.norm(candidate)+B*np.linalg.norm(ell-a@candidate)+scale*s[0]+float(c.ravel()@np.hypot(candidate[:2],candidate[2:]))
        assert fit['model_objective']<=exact.value+1e-6
        assert objective>=exact.value-1e-6
        assert fit['support_multipliers'].sum()<=scale*(1+1e-12)
        # Exported full-space support obeys the true spectral penalty globally.
        other=rng.normal(size=m);fo=np.einsum('j,jab->ab',other,fields)
        assert fit['pose_support']@other<=scale*np.linalg.norm(fo,2)+1e-10
        assert fit['cubic_support']@other<=float(c.ravel()@np.hypot(other[:2],other[2:]))+1e-10
        history.append(objective-fit['model_objective']);cuts.extend([support,-support])
        if history[-1]<1e-6:break
    assert history[-1]<1e-5
    np.testing.assert_allclose(objective,exact.value,rtol=1e-5,atol=1e-6)


def test_density_gram_rejects_material_indefiniteness():
    import pytest
    with pytest.raises(ValueError):psd_root(np.diag([1.,-.1]))


def test_cone_dual_uses_interior_subgradient_at_zero_weight_pair():
    ell=np.array([2.,.1,0.,0.]);joint=np.column_stack([ell,-np.eye(4)])
    c=np.array([[.05,.8]]);z=.1;B=.5
    fit=solve_pose_cut_problem(np.eye(4),joint.T@joint,c,np.zeros((0,4)),z,B,.1)
    w=fit['weights'];gr=fit['cubic_support']
    assert np.hypot(w[1],w[3])<1e-6
    assert np.hypot(gr[1],gr[3])<.7  # Strictly inside the 0.8 dual ball.
    gradient=z*w/np.linalg.norm(w)+B*(w-ell)/np.linalg.norm(w-ell)+gr
    assert np.linalg.norm(gradient)<1e-4
