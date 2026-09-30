import numpy as np
from scipy.optimize import minimize
from fourier_splats.uq_mixture_curvature import (
    CurvatureGaussianOrbit,euler_second_derivatives,quadratic_box_upper)
from fourier_splats.uq_continuous_mixture import euler_rotation_jacobian


def fixture():
    rng=np.random.default_rng(9871)
    orbit=CurvatureGaussianOrbit(np.pad(rng.normal(size=(8,2)),((0,0),(0,1))),
        rng.normal(size=(5,3)),np.linspace(.5,1.4,5),rng.normal(size=10),
        rng.uniform(-1,1,(3,8)),rng.normal(size=(3,16)))
    return rng,orbit


def test_second_rotation_and_log_derivatives():
    rng,orbit=fixture()
    for _ in range(5):
        angle=rng.normal(size=3);dd=euler_second_derivatives(angle)
        _,grad,hess,_,_=orbit.log_kernel_derivatives(angle)
        for j in range(3):
            step=np.eye(3)[j]*1e-5
            dp=euler_rotation_jacobian(angle+step)[1]
            dm=euler_rotation_jacobian(angle-step)[1]
            np.testing.assert_allclose((dp-dm)/2e-5,dd[:,j],atol=2e-10)
            fp,gp,*_=orbit.log_kernel_derivatives(angle+step)
            fm,gm,*_=orbit.log_kernel_derivatives(angle-step)
            np.testing.assert_allclose((fp-fm)/2e-5,grad[:,j],atol=2e-8)
            np.testing.assert_allclose((gp-gm)/2e-5,hess[:,:,j],atol=2e-8)


def test_quadratic_majorant_including_early_stopping():
    rng=np.random.default_rng(7413)
    for _ in range(16):
        g=rng.normal(size=3);a=rng.normal(size=(3,3));h=a+a.T
        box=rng.uniform(.01,1,3)
        bound=quadratic_box_upper(g,h,box)
        unfinished=quadratic_box_upper(g,h,box,sweeps=0)
        points=rng.uniform(-1,1,(1000,3))*box
        values=points@g+.5*np.einsum('ni,ij,nj->n',points,h,points)
        assert max(values)<=bound+1e-10 and bound<=unfinished+1e-10
        for p in points[:5]:
            opt=minimize(lambda x:-(g@x+.5*x@h@x),p,jac=lambda x:-(g+h@x),
                bounds=list(zip(-box,box)),method='L-BFGS-B')
            assert -opt.fun<=bound+1e-8
    diagonal=np.array([.4,1.1,3.]);g=np.array([.1,2.,-1.]);box=np.array([.7,.2,.9])
    point=np.clip(g/diagonal,-box,box)
    exact=g@point-.5*np.sum(diagonal*point**2)
    np.testing.assert_allclose(quadratic_box_upper(g,-np.diag(diagonal),box),exact,atol=1e-12)


def test_nonlinear_curvature_envelopes_and_third_remainder():
    rng,orbit=fixture()
    for diameter in [.001,.03,.2,1.]:
        center=rng.uniform(-3,3,3);half=rng.uniform(.2,1,3)*diameter
        cell=orbit.cell(center-half,center+half)
        ell,g,h,*_=orbit.log_kernel_derivatives(center)
        assert np.all(cell['envelope']<=cell['first_order_envelope'])
        for _ in range(100):
            d=rng.uniform(-1,1,3)*half
            actual=orbit.log_kernel(center+d)
            approximation=ell+g@d+.5*np.einsum('j,njk,k->n',d,h,d)
            assert np.all(np.abs(actual-approximation)<=cell['third_log_remainder']+1e-9)
            assert np.all(actual<=cell['quadratic_envelope']+1e-9)
            assert np.all(actual<=cell['envelope']+1e-9)
