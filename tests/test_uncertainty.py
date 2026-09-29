import numpy as np
import pytest
from scipy.stats import norm
from fourier_splats.uncertainty import optimize_certificate, bias_aware_half_width,compress_nuisance_group


def test_folded_normal_critical_values():
    for sd in [0.01, 1, 50]:
        for bias in [0, 0.2, 3, 100]:
            q = bias_aware_half_width(sd, bias)
            coverage = norm.cdf((q-bias)/sd)-norm.cdf((-q-bias)/sd)
            assert abs(coverage-0.95) < 1e-10
    assert bias_aware_half_width(0, 3) == 3


def test_scalar_exact_optimum_and_unidentified_target():
    # min z|w|+B|1-a*w| has an analytic solution at w=0 or w=1/a.
    for a in [0, 0.1, 10]:
        fit = optimize_certificate(np.array([[a]]), np.zeros((1,1,0)), [1], 2, maxiter=400, rtol=1e-5)
        expected = min(2, norm.ppf(0.975)/a) if a else 2
        assert fit.objective == pytest.approx(expected, rel=1e-5, abs=1e-7)
        assert fit.dual_lower_bound <= expected + 1e-8
    assert fit.noise_sd > 0


def test_primal_and_dual_against_independent_conic_solver():
    cp = pytest.importorskip("cvxpy")
    rng = np.random.default_rng(19)
    n,m,p,q = 8,6,9,2
    a = rng.normal(size=(n*m,p))
    j = rng.normal(size=(n,m,q))
    ell = rng.normal(size=p)
    B,eta,gamma = 3.0,0.8,0.12
    w = cp.Variable(n*m)
    cost = norm.ppf(0.975)*cp.norm(w) + B*cp.norm(ell-a.T@w)
    for i in range(n):
        wi=w[i*m:(i+1)*m]
        cost += eta*cp.norm(j[i].T@wi) + gamma*cp.norm(wi)
    oracle=cp.Problem(cp.Minimize(cost))
    optimum=oracle.solve(solver="CLARABEL",tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10)
    fit=optimize_certificate(a,j,ell,B,eta,gamma,maxiter=400,rtol=1e-5,smoothing=1e-9)
    assert fit.dual_lower_bound <= optimum + 1e-7
    assert fit.objective >= optimum - 1e-7
    assert fit.objective-optimum < 1e-4*optimum
    assert fit.converged


def test_coverage_for_adversarial_bias_and_nuisance():
    rng=np.random.default_rng(29)
    n,m,p,q=10,8,7,2
    a=rng.normal(size=(n*m,p))
    j=rng.normal(size=(n,m,q))
    ell=rng.normal(size=p)
    B,eta,gamma=2.0,0.6,0.2
    fit=optimize_certificate(a,j,ell,B,eta,gamma,maxiter=200)
    w=fit.weights
    delta=a.T@w-ell
    delta*=B/max(np.linalg.norm(delta),1e-30)
    jt=np.einsum('nmq,nm->nq',j,w.reshape(n,m))
    u=jt*eta/np.maximum(np.linalg.norm(jt,axis=1,keepdims=True),1e-30)
    r=w.reshape(n,m)*gamma/np.maximum(np.linalg.norm(w.reshape(n,m),axis=1,keepdims=True),1e-30)
    mean=a@delta+(np.einsum('nmq,nq->nm',j,u)+r).ravel()
    achieved_bias=w@mean-ell@delta
    assert achieved_bias == pytest.approx(fit.bias,rel=1e-9,abs=1e-9)
    # Simulate the exact scalar marginal; no pseudo-independent voxel trials.
    error=achieved_bias+rng.normal(scale=fit.noise_sd,size=100000)
    coverage=np.mean(np.abs(error)<=fit.half_width)
    assert abs(coverage-0.95)<0.002


def test_sparse_and_dense_agree():
    from scipy.sparse import csr_matrix
    rng=np.random.default_rng(31)
    a=rng.normal(size=(30,5)); ell=rng.normal(size=5)
    j=rng.normal(size=(5,6,1))
    dense=optimize_certificate(a,j,ell,2,0.3,0.1,maxiter=200)
    sparse=optimize_certificate(csr_matrix(a),j,ell,2,0.3,0.1,maxiter=200)
    assert sparse.objective == pytest.approx(dense.objective,rel=1e-8)
    assert np.allclose(sparse.weights,dense.weights,atol=1e-7)


@pytest.mark.parametrize('extra_rank',[1,9])
def test_multiple_groups_against_conic_solver(extra_rank):
    cp=pytest.importorskip('cvxpy')
    rng=np.random.default_rng(104)
    n,m,p,q=5,7,4,2
    a=rng.normal(size=(n*m,p));j=rng.normal(size=(n,m,q))
    k=rng.normal(size=(n,m,extra_rank));ell=rng.normal(size=p)
    B,eta,radius,gamma=4,.3,.12,.04
    w=cp.Variable(n*m)
    cost=norm.ppf(.975)*cp.norm(w)+B*cp.norm(ell-a.T@w)
    for i in range(n):
        wi=w[i*m:(i+1)*m]
        cost+=eta*cp.norm(j[i].T@wi)+radius*cp.norm(k[i].T@wi)+gamma*cp.norm(wi)
    optimum=cp.Problem(cp.Minimize(cost)).solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10)
    fit=optimize_certificate(a,j,ell,B,eta,gamma,extra_nuisance_groups=[(k,radius)],maxiter=500,rtol=2e-5)
    assert fit.dual_lower_bound<=optimum+1e-7
    assert abs(fit.objective-optimum)<1e-4*optimum
    assert fit.nuisance_bias==pytest.approx(sum(fit.nuisance_group_biases))


def test_wide_group_gram_reduction():
    rng=np.random.default_rng(418)
    d=rng.normal(size=(4,9,170));factor=compress_nuisance_group(d)
    w=rng.normal(size=(4,9))
    original=np.linalg.norm(np.einsum('nmp,nm->np',d,w),axis=1)
    compressed=np.linalg.norm(np.einsum('nmp,nm->np',factor,w),axis=1)
    assert np.all(compressed>=original-1e-12)
    assert np.allclose(compressed,original,rtol=1e-10)
