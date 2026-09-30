import numpy as np
import cvxpy as cp
import pytest
from scipy.sparse.linalg import LinearOperator
from fourier_splats.uq_two_pose_modulus import PhaseRotatedGram,two_pose_witness,optimize_two_pose_modulus
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram


class DenseGram:
    def __init__(self,A,ell):self.A=A;self.ell=ell;self.shape=(len(A),len(A));self.G=A@A.T
    def matvec(self,w):return self.G@w
    def target(self,*_):return self.A@self.ell,float(self.ell@self.ell)
    def quadrature_error(self,w):return {'squared_field_norm':0.,'gram_action_norm':0.}
    def preconditioner(self,ridge):
        return LinearOperator(self.shape,matvec=lambda v:np.linalg.solve(self.G+ridge*np.eye(len(self.A)),v))


def test_two_pose_witness_is_in_original_class_and_has_small_mean_distance():
    rng=np.random.default_rng(610331);ell=rng.normal(size=7);pilot=.1*rng.normal(size=7)
    matrices=[rng.normal(size=(4,7)) for _ in range(2)];grams=[DenseGram(A,ell) for A in matrices]
    B=2.;tau=.4;means=[A@pilot for A in matrices];targets=[A@ell for A in matrices]
    for size in [.001,1.,1000.]:
        w=size*rng.normal(size=4);record=two_pose_witness(grams,targets,ell@ell,means,w,B,np.linalg.norm(pilot),tau)
        a0,a1=record['residual_amplitudes'];scale=record['common_density_scale']
        rho0=scale*(pilot-a0*(ell-matrices[0].T@w));rho1=scale*(pilot+a1*(ell-matrices[1].T@w))
        assert max(np.linalg.norm(rho0-pilot),np.linalg.norm(rho1-pilot))<=B+1e-12
        assert np.linalg.norm(matrices[1]@rho1-matrices[0]@rho0)<=tau+1e-10
        np.testing.assert_allclose(ell@(rho1-rho0),record['signed_feature_gap'],atol=1e-10)


def test_two_pose_modulus_brackets_independent_conic_solution():
    rng=np.random.default_rng(610332);ell=rng.normal(size=6);pilot=.1*rng.normal(size=6)
    matrices=[rng.normal(size=(4,6)) for _ in range(2)];B=1.;tau=.7
    grams=[DenseGram(A,ell) for A in matrices];means=[A@pilot for A in matrices]
    initial=np.linalg.solve(sum(A@A.T for A in matrices)+np.eye(4),sum(A@ell for A in matrices))
    fit=optimize_two_pose_modulus(grams,None,None,None,means,initial,B,np.linalg.norm(pilot),tau,maxiter=100,rtol=1e-5)
    import json
    json.dumps({key:value for key,value in fit.items() if key not in ['weights','upper_weights']})
    rho0=cp.Variable(6);rho1=cp.Variable(6)
    problem=cp.Problem(cp.Maximize(ell@(rho1-rho0)),[cp.norm(rho0-pilot)<=B,cp.norm(rho1-pilot)<=B,
                 cp.norm(matrices[1]@rho1-matrices[0]@rho0)<=tau])
    optimum=problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10)
    assert fit['feasible_lower']<=optimum+1e-7<=fit['dual_upper']+1e-7
    assert fit['relative_modulus_gap']<1e-4


def test_phase_wrapper_matches_independent_real_fourier_matrix():
    rng=np.random.default_rng(610333);k=rng.normal(size=(2,3,3));ctf=rng.normal(size=(2,3));translations=rng.normal(size=(2,3))
    base=QuadratureObservationGram(k,ctf,.8,order=14,preconditioner_rank=0);op=PhaseRotatedGram(base,translations)
    xyz=np.stack(base.nodes,axis=1)
    phase=-2*np.pi*(np.einsum('nqa,xa->nqx',k,xyz)+translations[:,:,None])
    A=np.concatenate([np.cos(phase),np.sin(phase)],axis=1).reshape(12,-1)
    transfer=np.concatenate([ctf/.8,ctf/.8],axis=1).ravel()
    A*=transfer[:,None]*np.sqrt(base.weights)[None,:]
    w=rng.normal(size=12)
    np.testing.assert_allclose(op.matvec(w),A@(A.T@w),atol=1e-10,rtol=1e-10)
    np.testing.assert_allclose(op.rotate(op.rotate(w),True),w,atol=1e-14)


def test_inconsistent_residual_and_nonfinite_inputs_raise():
    g=DenseGram(np.eye(2),np.ones(2));w=np.ones(2)
    with pytest.raises(FloatingPointError,match='residual square'):
        two_pose_witness([g,g],[w,w],0.,[np.zeros(2)]*2,w,1.,0.,.7)
    with pytest.raises(FloatingPointError,match='inputs'):
        two_pose_witness([g,g],[w,w],2.,[np.zeros(2)]*2,w*np.nan,1.,0.,.7)


def test_crude_quadrature_pads_cover_exact_continuous_witness_and_phase_paths():
    from fourier_splats.uq_continuous import ContinuousObservationGram,cell_forward
    rng=np.random.default_rng(610334);n=2;nq=3;ctf=rng.normal(size=(n,nq));noise=.8
    pilot=rng.normal(size=4**3);pilot*=.3/np.linalg.norm(pilot);w=.1*rng.normal(size=2*n*nq)
    approximate=[];exact=[];means=[];targets=[]
    for _ in range(2):
        k=.7*rng.normal(size=(n,nq,3));t=.1*rng.normal(size=(n,nq))
        base=QuadratureObservationGram(k,ctf,noise,order=3,preconditioner_rank=0)
        op=PhaseRotatedGram(base,t);precise=PhaseRotatedGram(ContinuousObservationGram(k,ctf,noise),t)
        approximate.append(op);exact.append(precise);a,norm2=op.target([[0,0,0]],[1],.2);targets.append(a)
        means.append(op.rotate(cell_forward(k,ctf,pilot,4,noise)))
        # Explicit real rotation, including the packed per-particle row order.
        R=np.zeros((2*n*nq,2*n*nq))
        for i in range(n):
            for q in range(nq):
                j=i*2*nq+q;h=j+nq;c=np.cos(2*np.pi*t[i,q]);s=np.sin(2*np.pi*t[i,q])
                R[j,j]=R[h,h]=c;R[j,h]=s;R[h,j]=-s
        original,_=base.target([[0,0,0]],[1],.2)
        np.testing.assert_allclose(a,R@original,atol=1e-12)
        np.testing.assert_allclose(op.preconditioner(1.2)@w,R@(base.preconditioner(1.2)@(R.T@w)),atol=1e-12)
        e=op.quadrature_error(w)
        assert e['squared_field_norm']>0 and e['gram_action_norm']>0
        assert abs(w@(op.matvec(w)-precise.matvec(w)))<=e['squared_field_norm']+1e-12
        assert np.linalg.norm(op.matvec(w)-precise.matvec(w))<=e['gram_action_norm']+1e-12
    record=two_pose_witness(approximate,targets,norm2,means,w,1.,.3,.7)
    for j in range(2):
        actual2=norm2-2*w@targets[j]+w@exact[j].matvec(w)
        assert actual2<=record['residual_norms_upper'][j]**2+1e-12
    exact_mean=means[1]-means[0]+sum(a*(t-g.matvec(w)) for a,t,g in zip(record['residual_amplitudes'],targets,exact))
    assert np.linalg.norm(exact_mean)<=record['unscaled_mean_distance_upper']+1e-12
    assert record['common_density_scale']*np.linalg.norm(exact_mean)<=record['mean_distance_upper']+1e-12


def test_wrapped_pilot_mean_matches_nonlinear_constant_cell_projection():
    from fourier_splats.uq_continuous import cell_forward
    from fourier_splats.uq_continuous_pose import pose_cell_forward,perturbed_geometry
    rng=np.random.default_rng(610335);k=rng.normal(size=(2,3,3));q=rng.normal(size=(2,3,2));ctf=rng.normal(size=(2,3))
    pilot=rng.normal(size=4**3);pose=rng.normal(size=(2,5));pose/=np.linalg.norm(pose,axis=1)[:,None]
    rotated,t=perturbed_geometry(k,q,pose,.04,.007)
    op=PhaseRotatedGram(QuadratureObservationGram(rotated,ctf,.8,order=4,preconditioner_rank=0),t)
    np.testing.assert_allclose(op.rotate(cell_forward(rotated,ctf,pilot,4,.8)),
                              pose_cell_forward(k,q,ctf,pilot,4,.8,pose,.04,.007),atol=1e-12)


def test_separate_upper_and_lower_witness_weights_are_preserved():
    rng=np.random.default_rng(29);ell=rng.normal(size=6);pilot=.1*rng.normal(size=6)
    matrices=[rng.normal(size=(4,6)) for _ in range(2)];grams=[DenseGram(A,ell) for A in matrices]
    initial=3*rng.normal(size=4);means=[A@pilot for A in matrices];targets=[A@ell for A in matrices]
    fit=optimize_two_pose_modulus(grams,None,None,None,means,initial,1.,np.linalg.norm(pilot),.7,maxiter=6,rtol=1e-10,cg_maxiter=1)
    assert fit['best_lower_iteration']!=fit['best_upper_iteration']
    upper=two_pose_witness(grams,targets,ell@ell,means,fit['upper_weights'],1.,np.linalg.norm(pilot),.7)
    lower=two_pose_witness(grams,targets,ell@ell,means,fit['weights'],1.,np.linalg.norm(pilot),.7)
    np.testing.assert_allclose(upper['fixed_pair_modulus_upper'],fit['dual_upper'])
    np.testing.assert_allclose(lower['feature_gap'],fit['feasible_lower'])
