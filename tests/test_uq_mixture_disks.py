import numpy as np
import pytest
from fourier_splats.uq_mixture_disks import (
    DiskGaussianOrbit,disk_residual_dual,validate_independent_plane)


def test_disk_duals_against_independent_conic_program():
    import cvxpy as cp
    rng=np.random.default_rng(941)
    r=rng.normal(size=(4,12));j=rng.normal(size=(4,12,3))
    eps=rng.uniform(.01,.4,size=(4,6));h=np.array([.12,.2,.31])
    fit=disk_residual_dual(r,j,eps,h,sweeps=400)
    for i in range(4):
        v=cp.Variable(3);e=cp.Variable((6,2))
        packed=cp.hstack([e[:,0],e[:,1]])
        problem=cp.Problem(cp.Minimize(cp.sum_squares(r[i]-j[i]@v-packed)),
            [v<=h,v>=-h,cp.norm(e,axis=1)<=eps[i]])
        optimum=problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,
            tol_gap_rel=1e-10,tol_feas=1e-10,max_iter=500)
        assert problem.status=='optimal'
        assert fit['lower_squared'][i]<=optimum+1e-7
        assert fit['primal_squared'][i]>=optimum-1e-7
        assert fit['gap'][i]<2e-6
        u=fit['dual_vectors'][i]
        replay=(2*r[i]@u-u@u-2*eps[i]@np.hypot(u[:6],u[6:])
                -2*h@np.abs(j[i].T@u))
        np.testing.assert_allclose(replay,fit['lower_squared'][i],atol=1e-12)


def test_disks_preserve_frequency_error_allocation():
    # An image-wide radius spends error on a coordinate whose Taylor error is
    # actually zero. Per-frequency disks retain that observed residual.
    r=np.array([[3.,4.,0.,0.]]);j=np.zeros((1,4,3))
    fit=disk_residual_dual(r,j,np.array([[3.,0.]]),np.ones(3),sweeps=0)
    np.testing.assert_allclose(fit['lower_squared'],[16.])
    np.testing.assert_allclose(fit['primal_squared'],[16.])
    assert (np.linalg.norm(r)-3)**2==4.
    zero=disk_residual_dual(r,j,np.zeros((1,2)),np.zeros(3),sweeps=0)
    np.testing.assert_allclose(zero['lower_squared'],[25.])


def test_frequency_guards_and_nonlinear_cell_witnesses():
    for invalid in [np.array([[0.,0.,0.]]),np.array([[1.,2.,0.],[1.,2.,0.]]),
                    np.array([[1.,2.,0.],[-1.,-2.,0.]])]:
        with pytest.raises(ValueError):validate_independent_plane(invalid)
    rng=np.random.default_rng(942)
    plane=np.column_stack([rng.normal(size=(8,2)),np.zeros(8)])
    orbit=DiskGaussianOrbit(plane,rng.normal(size=(4,3)),.7,
        rng.normal(size=8),rng.normal(size=(3,8)),rng.normal(size=(3,16)))
    for half in [.0001,.02,.1,.6]:
        middle=np.array([1.1,.8,2.3]);lo=middle-half;hi=middle+half
        cell=orbit.cell(lo,hi)
        assert np.all(cell['envelope']<=cell['previous_envelope'])
        for angle in rng.uniform(lo,hi,size=(80,3)):
            assert np.all(orbit.log_kernel(angle)<=cell['envelope']+1e-9)
    point=orbit.cell(np.ones(3),np.ones(3))
    np.testing.assert_allclose(point['envelope'],point['exact'],atol=1e-10)
