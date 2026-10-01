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


def test_singular_galerkin_cache_retains_valid_old_bound(monkeypatch):
    g=np.diag([2.,1.]);b=np.array([.1,1.]);q=np.array([[1.,1.],[0.,0.]])
    result=quadratic_dual_upper(b,1.,2.,q,g@q)
    assert result['uses_old_cross_fallback'] and result['failed_evaluations']
    assert result['quadratic_upper']==result['old_cross_quadratic_upper']
    # LAPACK may numerically solve some rank-deficient shifts. Force the
    # actual all-solve-fail branch instead of assuming each factorization fails.
    def fail(*args,**kwargs):raise np.linalg.LinAlgError('injected solve failure')
    monkeypatch.setattr(np.linalg,'solve',fail)
    result=quadratic_dual_upper(b,1.,2.,np.eye(2),g)
    assert result['selected'] is None and len(result['failed_evaluations'])==61
    assert result['quadratic_upper']==result['old_cross_quadratic_upper']


def test_near_hard_log_root_is_accurate_before_boundary_extension():
    for epsilon in [1e-8,1e-10,1e-12]:
        b=np.array([epsilon,1.]);g=np.diag([5.,2.]);L=2.
        v,r=dense_quadratic_witness(g,b,L)
        assert not r['hard_case']
        np.testing.assert_allclose(r['pre_boundary_norm'],L,rtol=2e-12)
        raw=v*r['pre_boundary_norm']/np.linalg.norm(v)
        expected=20+1/3+2*np.sqrt(4-1/9)*epsilon
        np.testing.assert_allclose(raw@g@raw-2*b@raw,expected,atol=1e-10,rtol=0)


def test_truncated_cache_and_loose_spectral_upper_report_realistic_slack():
    rng=np.random.default_rng(953107);m=rng.normal(size=(23,23));g=m.T@m/23;b=rng.normal(size=23)
    q,gq,cache=krylov_cache(lambda x:g@x,b,steps=3,extras=[rng.normal(size=23)])
    U=1.5*np.linalg.eigvalsh(g)[-1]
    result=quadratic_dual_upper(b,1.4,U,q,gq);v,_=dense_quadratic_witness(g,b,1.4)
    assert result['quadratic_upper']>=v@g@v-2*b@v-1e-10
    assert result['quadratic_upper']<=result['old_cross_quadratic_upper']
    assert [row['origin'] for row in cache['accepted_chains']]==['cross','extra_0','cross']
    print('LOOSE_TRUNCATED_AUDIT',result['quadratic_upper'],result['old_cross_quadratic_upper'],result['uses_old_cross_fallback'])
