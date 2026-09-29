import numpy as np
from fourier_splats.uq_continuous_pose import pose_derivative_fields, derivative_coefficient_bounds
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_random_spectral import gaussian_power_upper


def test_matrix_free_pose_fields_and_adjoint_against_dense_derivatives():
    rng = np.random.default_rng(609581)
    n, nq = 3, 7
    k = rng.normal(size=(n,nq,3)); q = rng.normal(size=(n,nq,2))
    ctf = rng.normal(size=(n,nq)); w = rng.normal(size=2*n*nq)
    angle = np.array([.01,.03,.07]); shift = np.array([.001,.002,.004])
    op = PolynomialPoseFieldOperator(k,q,ctf,w,.7,angle,shift,order=6,backend='direct',block_particles=2)
    f = np.empty(op.shape)
    for i in range(n):
        a,b=pose_derivative_fields(k[i],q[i],op.c[i],op.xyz,angle[i],shift[i],backend='direct')
        f[:,5*i:5*(i+1)] = op.sqrt_quad[:,None]*a
        f[:,5*n+15*i:5*n+15*(i+1)] = .5*op.sqrt_quad[:,None]*b
    u=rng.normal(size=20*n);v=rng.normal(size=len(op.xyz))
    np.testing.assert_allclose(op.matvec(u),f@u,rtol=1e-12,atol=2e-14)
    np.testing.assert_allclose(op.rmatvec(v),f.T@v,rtol=1e-12,atol=2e-14)
    fast=PolynomialPoseFieldOperator(k,q,ctf,w,.7,angle,shift,order=6,backend='nufft',block_particles=2)
    np.testing.assert_allclose(fast.matvec(u),f@u,rtol=1e-9,atol=2e-11)
    np.testing.assert_allclose(fast.rmatvec(v),f.T@v,rtol=1e-9,atol=2e-11)
    scaling=op.establish_group_scaling()
    scaled=f/op.denominators[None]
    np.testing.assert_allclose(op.spatial_gram(v),scaled@(scaled.T@v),rtol=1e-12,atol=2e-13)
    np.testing.assert_allclose(scaling['quadrature_gram_trace'],np.sum(scaled**2),rtol=1e-12)
    amplitudes=np.hypot(w.reshape(n,2*nq)[:,:nq],w.reshape(n,2*nq)[:,nq:])
    masses=np.einsum('naj,nj->na',op.coefficient_bound_matrix(),amplitudes)
    expected=[]
    for i in range(n):
        m1,m2=derivative_coefficient_bounds(k[i],q[i],op.c[i],angle[i],shift[i])
        den=np.concatenate([op.denominators[5*i:5*(i+1)],op.denominators[5*n+15*i:5*n+15*(i+1)]])
        expected.append(np.concatenate([m1,.5*m2])/den)
    np.testing.assert_allclose(masses,expected,rtol=1e-12)
    gradient=op.weight_gradient(u,v)
    direction=rng.normal(size=w.shape);delta=1e-5
    op.set_weights(w+delta*direction); plus=v@op.matvec(u)
    op.set_weights(w-delta*direction); minus=v@op.matvec(u)
    np.testing.assert_allclose((plus-minus)/(2*delta),gradient@direction,rtol=1e-9,atol=1e-10)
    op.set_weights(w)
    envelope=op.establish_coefficient_scaling()
    scaled=f/op.denominators[None]
    assert envelope['quadrature_gram_trace_upper']>=np.sum(scaled**2)
    np.testing.assert_allclose(op.spatial_gram(v),scaled@(scaled.T@v),rtol=1e-12,atol=2e-13)


def test_random_power_bound_matches_direct_matrix_powers_and_stated_event():
    rng=np.random.default_rng(609582)
    q,_=np.linalg.qr(rng.normal(size=(12,12)))
    eigenvalues=np.linspace(.1,2,12); matrix=(q*eigenvalues)@q.T
    seed=17;probes=4;iterations=12;delta=1e-5
    result=gaussian_power_upper(lambda v:matrix@v,12,delta,probes,iterations,seed,trace_upper=np.trace(matrix))
    g=np.random.default_rng(seed).normal(size=(probes,12))
    event=np.max(abs(g@q[:,-1]))>=result['projection_threshold']
    assert event
    for p,row in enumerate(result['history'],start=1):
        direct=(np.max(np.linalg.norm(g@np.linalg.matrix_power(matrix,p),axis=1))/result['projection_threshold'])**(1/p)
        np.testing.assert_allclose(row['power_upper'],direct,rtol=1e-12)
        assert row['best_upper']>=2-1e-12
    assert result['rayleigh_lower']<=2+1e-12
    assert result['eigenvalue_upper']<2.8
