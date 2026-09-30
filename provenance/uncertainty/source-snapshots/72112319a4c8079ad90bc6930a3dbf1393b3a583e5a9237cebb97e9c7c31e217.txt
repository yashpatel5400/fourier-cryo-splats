import numpy as np
from scipy.optimize import lsq_linear
from scipy.spatial.transform import Rotation
from scipy.special import logsumexp
from fourier_splats.uq_continuous_mixture import (
    FourierGaussianOrbit, euler_rotation_jacobian, linear_box_residual_dual,
    continuous_mixture_upper)
from fourier_splats.uq_mixture_validation import mixture_envelope_upper


def fixture(seed=921):
    rng=np.random.default_rng(seed)
    plane=np.pad(rng.normal(size=(9,2)),((0,0),(0,1)))
    centers=rng.normal(size=(5,3))
    coefficient=rng.normal(size=10)
    transfer=rng.uniform(-1,1,(4,9))
    return rng,FourierGaussianOrbit(plane,centers,np.linspace(.5,1.4,5),
        coefficient,transfer,rng.normal(size=(4,18)))


def test_euler_derivatives_and_image_jacobians_independently():
    rng,orbit=fixture()
    for _ in range(8):
        angle=rng.uniform(-4,4,3);r,dr=euler_rotation_jacobian(angle)
        np.testing.assert_allclose(r,Rotation.from_euler('ZYZ',angle).as_matrix(),atol=1e-14)
        m,j,_=orbit.mean_jacobian(angle)
        for a in range(3):
            step=np.eye(3)[a]*1e-5
            rp,_=euler_rotation_jacobian(angle+step)
            rm,_=euler_rotation_jacobian(angle-step)
            np.testing.assert_allclose((rp-rm)/2e-5,dr[a],atol=3e-10)
            mp=orbit.mean_jacobian(angle+step)[0]
            mm=orbit.mean_jacobian(angle-step)[0]
            np.testing.assert_allclose((mp-mm)/2e-5,j[:,:,a],atol=2e-9)


def test_residual_dual_bounds_independent_bounded_least_squares():
    rng=np.random.default_rng(131)
    for scale in [.01,.3,3.]:
        for rank in [0,1,3]:
            j=rng.normal(size=(12,rank))@rng.normal(size=(rank,3))
            r=rng.normal(size=12);h=scale*rng.uniform(.1,1,3)
            exact=lsq_linear(j,r,bounds=(-h,h),tol=1e-13,lsq_solver='exact',max_iter=500)
            lower=linear_box_residual_dual(np.array(r@r),j.T@r,j.T@j,h,sweeps=32)
            assert lower>=0 and lower<=2*exact.cost+2e-10
            # Even an intentionally unfinished numerical solve remains a dual bound.
            unfinished=linear_box_residual_dual(np.array(r@r),j.T@r,j.T@j,h,sweeps=0)
            assert unfinished<=2*exact.cost+2e-10


def test_rotation_box_envelopes_and_taylor_remainders():
    rng,orbit=fixture(187)
    for half_scale in [.0001,.03,.3,2.]:
        middle=rng.uniform(0,3,3);half=half_scale*rng.uniform(.3,1,3)
        box=orbit.cell(middle-half,middle+half)
        m,j,_=orbit.mean_jacobian(middle)
        for _ in range(60):
            delta=rng.uniform(-1,1,3)*half
            value=orbit.log_kernel(middle+delta)
            assert np.all(value<=box['envelope']+2e-10)
            actual=orbit.mean_jacobian(middle+delta)[0]-m-np.einsum('ndj,j->nd',j,delta)
            assert np.all(np.linalg.norm(actual,axis=1)<=box['remainder']+2e-11)
        if half_scale==.0001:
            assert np.max(box['envelope']-box['exact'])<.1


def assert_cover(lower,upper):
    # Volume plus pointwise checks; bisection history gives the analytic cover.
    assert abs(np.prod(upper-lower,axis=1).sum()-4*np.pi**3)<1e-9
    samples=np.random.default_rng(387).uniform(size=(1000,3))*[2*np.pi,np.pi,2*np.pi]
    for p in samples:
        assert np.sum(np.all((p>=lower)&(p<=upper),axis=1))==1


def test_continuous_dual_replay_and_independent_grid():
    rng,orbit=fixture(991)
    fit=continuous_mixture_upper(orbit,initial_bins=(2,1,2),max_splits=32,
        refit_every=16,mixture_iterations=30,tolerance=1e-10,wall_seconds=60)
    assert not fit['converged'] and fit['status']=='split_limit'
    assert_cover(fit['leaf_lower'],fit['leaf_upper'])
    assert_cover(fit['best_cover_lower'],fit['best_cover_upper'])
    log_z=fit['best_anchor_log_z']
    replay=log_z.sum()+orbit.n*(np.max(logsumexp(
        fit['best_cover_log_envelopes']-log_z[None,:],axis=1))-np.log(orbit.n))+orbit.normalization
    np.testing.assert_allclose(replay,fit['log_likelihood_upper'],atol=1e-10)
    count=fit['best_anchor_support_count']
    np.testing.assert_allclose(logsumexp(fit['support_log_kernels'][:count].T+
        np.log(fit['best_anchor_weights'])[None,:],axis=1),log_z,atol=1e-12)
    angles=rng.uniform(size=(1500,3))*[2*np.pi,np.pi,2*np.pi]
    values=np.array([orbit.log_kernel(a) for a in angles]).T
    grid=mixture_envelope_upper(values,max_iterations=100,tolerance=.01)
    assert grid.envelope_primal+orbit.normalization<=fit['log_likelihood_upper']+1e-9
    assert fit['log_likelihood_lower']<=fit['log_likelihood_upper']+1e-9
    # No splitting and a zero map are a useful exact continuous control.
    orbit.coefficients[:]=0.;orbit.abs_coefficient[:]=0.
    zero=continuous_mixture_upper(orbit,initial_bins=(1,1,1),max_splits=0,tolerance=1e-8)
    exact=-.5*np.sum(orbit.observed**2)+orbit.normalization
    assert zero['converged']
    np.testing.assert_allclose([zero['log_likelihood_lower'],zero['log_likelihood_upper']],exact,atol=1e-10)
