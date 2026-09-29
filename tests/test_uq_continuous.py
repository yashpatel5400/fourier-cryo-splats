import numpy as np
from numpy.polynomial.legendre import leggauss
from fourier_splats.uq_continuous import (gaussian_target_integrals,fourier_field_norm2,continuous_residual_norm,
                                        cell_forward,cell_adjoint,cell_target_coefficients,CellObservationOperator,
                                        ContinuousObservationGram,continuous_certificate)


def test_continuous_formulas_against_independent_volume_quadrature():
    rng=np.random.default_rng(609365);k=rng.normal(size=(2,5,3));ctf=rng.normal(size=(2,5));w=rng.normal(size=20)
    centers=np.array([[.08,0,.03],[-.11,.05,-.08]]);signs=np.array([1,-.7]);sigma=.12;noise=.8
    nodes,weights=leggauss(44);nodes/=2;weights/=2
    z,y,x=np.meshgrid(nodes,nodes,nodes,indexing='ij');xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],-1)
    quad=(weights[:,None,None]*weights[None,:,None]*weights[None,None,:]).ravel()
    target=sum(s*np.exp(-np.sum((xyz-mu)**2,axis=1)/(2*sigma*sigma))/(2*np.pi*sigma*sigma)**1.5 for mu,s in zip(centers,signs))
    phase=np.exp(2j*np.pi*(xyz@k.reshape(-1,3).T));wb=w.reshape(2,10);c=((wb[:,:5]+1j*wb[:,5:])*ctf/noise).ravel()
    field=(phase@c).real;ft,tt=gaussian_target_integrals(k,centers,signs,sigma)
    np.testing.assert_allclose(ft.ravel(),(quad*target)@phase.conj(),atol=2e-12)
    np.testing.assert_allclose(tt,quad@target**2,rtol=1e-12)
    norm2,_=fourier_field_norm2(k,c,chunk=3);np.testing.assert_allclose(norm2,quad@field**2,rtol=1e-12)
    audit=continuous_residual_norm(k,ctf,w,noise,centers,signs,sigma)
    np.testing.assert_allclose(audit['residual_norm2_unpadded'],quad@(target-field)**2,rtol=1e-12)
    for box in [8,16]:
        coefficients=rng.normal(size=box**3);adjoint=cell_adjoint(k,ctf,w,box,noise)
        operator=CellObservationOperator(k,ctf,box,noise)
        np.testing.assert_allclose(operator.adjoint(w),adjoint,atol=1e-9)
        np.testing.assert_allclose(operator.forward(coefficients),cell_forward(k,ctf,coefficients,box,noise),atol=1e-9)
        np.testing.assert_allclose(w@cell_forward(k,ctf,coefficients,box,noise),adjoint@coefficients,atol=1e-9)
        residual=cell_target_coefficients(box,centers,signs,sigma)-adjoint
        assert np.linalg.norm(residual)<=audit['residual_norm']+1e-9
    # For nested cells the orthogonal projection captures increasing energy.
    p8=cell_target_coefficients(8,centers,signs,sigma)-cell_adjoint(k,ctf,w,8,noise)
    p16=cell_target_coefficients(16,centers,signs,sigma)-cell_adjoint(k,ctf,w,16,noise)
    assert np.linalg.norm(p8)<=np.linalg.norm(p16)+1e-9


def test_continuous_gram_and_dual_against_independent_conic_problem():
    import cvxpy as cp
    from scipy.stats import norm
    rng=np.random.default_rng(609367);k=rng.normal(size=(2,4,3));ctf=rng.normal(size=(2,4));noise=.9
    centers=[[0,0,0]];signs=[1];width=.18;B=.7
    gram=ContinuousObservationGram(k,ctf,noise);a,lnorm2=gram.target(centers,signs,width)
    dense=np.column_stack([gram.matvec(v) for v in np.eye(16)])
    w=rng.normal(size=16);audit=continuous_residual_norm(k,ctf,w,noise,centers,signs,width)
    np.testing.assert_allclose(w@dense@w,audit['field_norm2'],rtol=1e-12)
    np.testing.assert_allclose(w@a,audit['cross_term'],rtol=1e-12)
    # The joint Hilbert Gram of [ell, A*_1,...,A*_m] gives an independent finite
    # representation containing every direction relevant to this one objective.
    joint=np.block([[np.array([[lnorm2]]),a[None]],[a[:,None],dense]])
    vals,vecs=np.linalg.eigh(joint);assert vals.min()>-1e-10
    factors=vecs*np.sqrt(np.maximum(vals,0))[None];ell=factors[0];A=factors[1:]
    variable=cp.Variable(16);z=norm.ppf(.975)
    problem=cp.Problem(cp.Minimize(z*cp.norm(variable)+B*cp.norm(ell-A.T@variable)))
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
    result=continuous_certificate(gram,centers,signs,width,B,rtol=1e-5,maxiter=200)
    assert result['converged']
    assert result['dual_lower_bound']<=problem.value+1e-6
    assert problem.value<=result['objective']+1e-6
    np.testing.assert_allclose(result['objective'],problem.value,rtol=2e-5)
