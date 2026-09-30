"""Compare residual bounds to independent exact solves and SDP duals."""
import numpy as np
import pytest
from fourier_splats.uq_trust_region import dense_quadratic_witness, krylov_cache, resolvent_upper, quadratic_dual_upper


def test_hard_and_near_hard_cases_with_analytic_answers():
    g = np.diag([5.,2.]); b = np.array([0.,1.]); radius=2.
    v, rec = dense_quadratic_witness(g,b,radius)
    assert rec['hard_case']
    np.testing.assert_allclose(v@g@v-2*b@v,5*radius**2+1/3,rtol=1e-12)
    assert np.linalg.norm(v)<=radius*(1+1e-15)
    near=np.array([1e-8,1.]); vv,_=dense_quadratic_witness(g,near,radius)
    assert vv@g@vv-2*near@vv>=rec['value']-1e-8
    for zero in [np.zeros(2), np.array([1.,0.])]:
        point,record=dense_quadratic_witness(g,zero,radius)
        assert np.isfinite(point).all() and np.linalg.norm(point)<=radius*(1+1e-15)
    with pytest.raises(ValueError): dense_quadratic_witness(-np.eye(2),b,radius)


def test_explicit_residual_controls_truncated_and_nonorthogonal_cache():
    rng=np.random.default_rng(953101); x=rng.normal(size=(19,19)); g=x.T@x/19
    b=rng.normal(size=19); U=np.linalg.eigvalsh(g)[-1]*1.1; lam=U+1.
    q,gq,_=krylov_cache(lambda x:g@x,b,steps=3)
    # Perturbing the basis does not invalidate the any-x residual identity.
    transform=np.array([[1.,.1,.2],[0.,.8,.2],[0.,0.,1.3]])
    q=q@transform;gq=gq@transform
    got=resolvent_upper(b,U,lam,q,gq)
    exact=float(b@np.linalg.solve(lam*np.eye(19)-g,b))
    assert got['upper']>=exact-1e-12
    assert got['residual_correction']>0
    full=resolvent_upper(b,U,lam,np.eye(19),g)
    np.testing.assert_allclose(full['upper'],exact,rtol=1e-12)
    empty=resolvent_upper(b,U,lam,np.zeros((19,0)),np.zeros((19,0)))
    assert empty['upper']>=exact
    with pytest.raises(ValueError):resolvent_upper(b,U,U,q,gq)


def test_joint_dual_against_independent_semidefinite_program():
    import cvxpy as cp
    rng=np.random.default_rng(953102)
    for hard in [False,True]:
        a=rng.normal(size=(9,6));g=a.T@a/9
        b=rng.normal(size=6)
        if hard:
            d,u=np.linalg.eigh(g);b=u[:,0]*.01
        L=1.4; U=float(np.linalg.eigvalsh(g)[-1])*(1+1e-10)
        # Independent standard trust-region SDP dual: s >= b'(lambda I-G)^-1 b.
        lam=cp.Variable(nonneg=True);s=cp.Variable()
        block=cp.bmat([[lam*np.eye(6)-g,b[:,None]],[b[None],cp.reshape(s,(1,1),order='C')]])
        problem=cp.Problem(cp.Minimize(lam*L*L+s),[block>>0])
        expected=problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10)
        assert problem.status=='optimal'
        v,record=dense_quadratic_witness(g,b,L)
        np.testing.assert_allclose(record['value'],expected,rtol=3e-8,atol=3e-8)
        q,gq,cache=krylov_cache(lambda x:g@x,b,steps=6,extras=[rng.normal(size=6)])
        got=quadratic_dual_upper(b,L,U,q,gq)
        assert got['quadratic_upper']>=expected-1e-7
        assert got['quadratic_upper']<=expected+2e-6
        assert got['quadratic_upper']<=got['old_cross_quadratic_upper']
        assert cache['orthogonality_defect']<1e-12
