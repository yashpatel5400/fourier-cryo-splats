"""Reduced design checked against independently assembled dense conic norms."""
import numpy as np
from scipy.linalg import block_diag
from test_uq_cubic_design_review import two_particle_fixture
from fourier_splats.uq_cubic_subspace import ReducedCubicDesign, orthonormal_basis, optimize_reduced_cubic


def dense_field(obj, weights):
    obj.op.set_weights(weights)
    return np.column_stack([obj.op.matvec(e) for e in np.eye(obj.op.shape[1])])


def reduced_fixture():
    obj, w = two_particle_fixture(); obj.smoothing = 0.
    rng = np.random.default_rng(710101)
    Q, record = orthonormal_basis(np.column_stack([w,rng.normal(size=w.size),rng.normal(size=w.size)]))
    assert record['rank']==3
    return obj,w,Q


def test_reduced_majorant_and_spectral_cuts_at_independent_weights():
    obj,w,Q=reduced_fixture(); reduced=ReducedCubicDesign(obj,Q)
    rng=np.random.default_rng(710102)
    x=Q.T@w
    sigma,cuts,_,_,_=reduced.spectral_directions(x,initial=np.ones(obj.op.shape[1]),modes=2,dense=True)
    assert abs(max(c@x for c in cuts)-sigma)<1e-10
    for _ in range(8):
        other=rng.normal(size=3); field=dense_field(obj,Q@other)
        u,s,_=np.linalg.svd(field,full_matrices=False)
        exact,_,_,_=obj.evaluate(Q@other,u[:,:1])
        majorant,_=reduced.value(other,s[0])
        assert majorant>=exact-1e-9
        assert max(abs(c@other) for c in cuts)<=s[0]+1e-10
    zero=reduced.spectral_directions(np.zeros(3),initial=np.ones(obj.op.shape[1]),modes=2,dense=True)
    assert zero[0]==0 and np.isfinite(zero[2]).all()


def test_pose_column_gram_oracle_agrees_with_dense_singular_value():
    obj,w,Q=reduced_fixture(); reduced=ReducedCubicDesign(obj,Q)
    start=np.random.default_rng(710105).normal(size=obj.op.shape[1])
    found=reduced.spectral_directions(Q.T@w,initial=start,modes=2,tolerance=1e-10,subspace=17)
    expected=np.linalg.svd(dense_field(obj,w),compute_uv=False)[0]
    np.testing.assert_allclose(found[0],expected,rtol=2e-10,atol=1e-12)


def independent_reduced_conic(obj,Q):
    import cvxpy as cp
    # Independently form all spatial matrices; avoid the production moment,
    # coefficient-gradient and spectral-oracle paths in this comparison.
    matrices=[dense_field(obj,q) for q in Q.T]
    left,sv,_=np.linalg.svd(np.concatenate(matrices,axis=1),full_matrices=False)
    right,sr,_=np.linalg.svd(np.concatenate([a.T for a in matrices],axis=1),full_matrices=False)
    left=left[:,sv>sv[0]*1e-13];right=right[:,sr>sr[0]*1e-13]
    fields=[left.T@m@right for m in matrices]
    def root(g):
        v,u=np.linalg.eigh(.5*(g+g.T));assert v.min()>-1e-9*max(1.,v.max())
        return np.sqrt(np.maximum(v,0.))[:,None]*u.T
    # The density Gram is exact analytic sinc integration, independent of NUFFT.
    k=obj.op.k.reshape(-1,3);t=obj.op.transfer.ravel()
    minus=np.prod(np.sinc(k[:,None]-k[None,:]),axis=-1)
    plus=np.prod(np.sinc(k[:,None]+k[None,:]),axis=-1)
    rr=.5*t[:,None]*t[None,:]*(minus+plus)
    ii=.5*t[:,None]*t[None,:]*(minus-plus)
    n,nq=obj.n,obj.nq
    # Source weights pack real/imaginary coordinates inside each particle.
    order=np.arange(2*n*nq).reshape(2,n,nq).transpose(1,0,2).ravel()
    gram=block_diag(rr,ii)[np.ix_(order,order)]+obj.eta*np.eye(2*n*nq)
    a=Q.T@obj.a
    aug=np.block([[Q.T@gram@Q,-a[:,None]],[-a[None],np.array([[obj.target_norm2]])]])
    q=Q.reshape(n,2,nq,3);amp=np.hypot(q[:,0],q[:,1])
    mass=np.einsum('naj,njp->nap',obj.coefficients,amp).reshape(-1,3)*np.sqrt(obj.kernel_error)
    pilot=[]
    for weight in Q.T:
        obj.op.set_weights(weight);pilot.append(obj.op.pair_moments(obj.moments))
    pilot=np.column_stack(pilot)
    x=cp.Variable(3);absolute=cp.Variable(3,nonneg=True)
    field=sum(x[j]*fields[j] for j in range(3))
    cost=obj.z*cp.norm(Q@x)+obj.B*cp.norm(root(aug)@cp.hstack([x,1.]))
    cost+=obj.B*obj.L*cp.norm(cp.hstack([cp.norm(field,2),mass@absolute]))+obj.L*cp.norm(pilot@x)
    for i in range(n):
        for d in range(obj.remainder.order):
            g=sum(q[i,c].T@obj.remainder.blocks[i,d,c]@q[i,c] for c in range(2))
            cost+=(obj.B+obj.P)*obj.remainder.multipliers[i,d]*cp.norm(root(g)@x)
    coefficient=(obj.B+obj.P)*np.einsum('nj,njp->p',obj.remainder.residual,amp)
    cost+=coefficient@absolute
    problem=cp.Problem(cp.Minimize(cost),[absolute>=x,absolute>=-x])
    value=problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9)
    assert problem.status=='optimal'
    return value


def test_cutting_design_matches_independent_dense_sdp():
    obj,w,Q=reduced_fixture();expected=independent_reduced_conic(obj,Q)
    result=optimize_reduced_cubic(obj,Q,w,optimization_seed=710103,certificate_seed=710104,
        max_evaluations=60,separation_tolerance=1e-7,modes=2,dense_oracle=True,power_iterations=50)
    assert abs(result['selected_approximate_objective']-expected)<3e-5
    assert result['sum_objective_upper']>=result['dual_lower_bound']
    assert result['spectral_upper_certificate']['seed']==710104
    assert result['stopped_on_restricted_guide_gap']
    print('SUBSPACE_CHECK',expected,result['selected_approximate_objective'],len(result['optimization_history']))
