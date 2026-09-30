"""Independent robust-norm SDP and direct continuous integrations."""
import numpy as np
from test_uq_cubic_design_review import two_particle_fixture
from fourier_splats.uq_cubic_subspace import orthonormal_basis
from fourier_splats.uq_joint_cubic_design import ReducedJointCubicDesign, design_joint_cubic, audit_joint_cubic


def setup():
    obj,w=two_particle_fixture();obj.L*=.03
    obj.remainder.multipliers*=.03;obj.remainder.residual*=.03
    rng=np.random.default_rng(953103);q,_=orthonormal_basis(np.column_stack([w,rng.normal(size=w.size)]))
    return obj,w,q


def independently_assemble_robust_sdp(red):
    import cvxpy as cp
    obj=red.obj;op=obj.op;fields=[]
    for q in red.Q.T:
        op.set_weights(q)
        fields.append(np.column_stack([op.matvec(e) for e in np.eye(op.shape[1])]))
    # Exact two-sided subspace compression of all small spatial/pose fields.
    left,s,_=np.linalg.svd(np.column_stack([red.target,red.nominal,*fields]),full_matrices=False)
    left=left[:,s>s[0]*1e-12]
    right,s,_=np.linalg.svd(np.concatenate([m.T for m in fields],axis=1),full_matrices=False)
    right=right[:,s>s[0]*1e-12]
    fields=[left.T@m@right for m in fields];a=left.T@red.nominal;ell=left.T@red.target
    x=cp.Variable(red.p);t=cp.Variable(nonneg=True);lam=cp.Variable(nonneg=True);ab=cp.Variable(red.p,nonneg=True)
    f=obj.L*sum(x[j]*fields[j] for j in range(red.p));h=ell-a@x
    n,k=fields[0].shape
    # Independent S-lemma epigraph, rescaled exactly to a unit pose ball.
    # This avoids the badly scaled lambda needed with a small physical L.
    block=cp.bmat([[cp.reshape(t-lam,(1,1),order='C'),np.zeros((1,k)),cp.reshape(h,(1,n),order='C')],
        [np.zeros((k,1)),lam*np.eye(k),-f.T],
        [cp.reshape(h,(n,1),order='C'),-f,t*np.eye(n)]])
    r=red.reduced
    cost=obj.z*cp.norm(red.Q@x)+obj.B*t+obj.L*cp.norm(r.pilot_root@x)
    cost+=sum(cp.norm(root@x) for root in r.remainder_roots)+r.residual_coefficients@ab
    problem=cp.Problem(cp.Minimize(cost),[ab>=x,ab>=-x,block>>0])
    val=problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9,max_iter=500)
    assert problem.status=='optimal'
    return val,x.value


def test_joint_cuts_and_dense_separation_agree_with_robust_sdp():
    obj,w,q=setup();args=([[.08,-.06,0.]],[1.],.24)
    red=ReducedJointCubicDesign(obj,q,*args);expected,solution=independently_assemble_robust_sdp(red)
    result=design_joint_cubic(obj,q,w,*args,optimization_seed=953104,evaluations=40,krylov_steps=20,
        design_seconds=120,separation_tolerance=1e-7,dense_oracle=True)
    assert abs(result['selected_guide']-expected)<1e-4
    assert result['selected_guide']<result['history'][0]['guide_objective']
    print('JOINT_SDP',expected,result['selected_guide'],len(result['history']))


def test_continuous_joint_audit_covers_direct_sampled_worst_case():
    obj,w,q=setup();args=([[.08,-.06,0.]],[1.],.24)
    red=ReducedJointCubicDesign(obj,q,*args)
    v,oracle=red.witness(q.T@w,start=np.ones(obj.op.shape[1]),dense=True)
    audited=audit_joint_cubic(obj,w,*args,certificate_seed=953105,power_iterations=50,krylov_steps=30)
    assert audited['bias_upper']<=audited['old_joint']['bias_upper']
    assert audited['bias_upper']>=obj.B*oracle['joint_norm_guide']-1e-7
    assert audited['spectral_upper_certificate']['seed']==953105
    assert 'relative_sum_gap' not in audited['old_triangle_audit']
    assert audited['resolvent']['quadratic_upper']<=audited['resolvent']['old_cross_quadratic_upper']
