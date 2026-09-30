import numpy as np
from numpy.polynomial.legendre import leggauss
from fourier_splats.uq_cell_moments import cell_fourier_moments,normalized_cell_moments,direct_cell_fourier_moments,check_cell_pose_pairings
from fourier_splats.uq_joint_bias import MONOMIALS,pair_pose_fourier_moments,joint_density_pose_bias
from fourier_splats.uq_continuous import cell_forward
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator


def test_local_cell_moments_at_zero_and_transition_against_separable_quadrature():
    box=24;k=box*np.array([0.,1e-8,-1e-5,.029/np.pi,.031/np.pi,.21,-.21])
    x,w=leggauss(30);x/=2*box;w/=2
    exact=np.stack([np.exp(2j*np.pi*k[:,None]*x)@(w*x**j) for j in range(3)],axis=-1)
    np.testing.assert_allclose(normalized_cell_moments(k,box),exact,rtol=1e-9,atol=2e-15)


def test_full_cell_moments_against_independent_composite_cell_integration():
    rng=np.random.default_rng(610291);box=4;coef=rng.normal(size=box**3)
    k=5*rng.normal(size=(2,3,3));k*=5/np.maximum(5,np.linalg.norm(k,axis=-1))[...,None]
    # A high-order rule inside EVERY cell integrates the discontinuous pilot;
    # sampling it on one global quadrature grid would not be an exact check.
    x,qw=leggauss(12);x/=2*box;qw/=2*box
    z,y,xx=np.meshgrid(x,x,x,indexing='ij');offset=np.stack([xx,y,z],axis=-1).reshape(-1,3)
    weights=(qw[:,None,None]*qw[None,:,None]*qw[None,None,:]).ravel()
    expected=np.zeros((6,10),complex);frequency=k.reshape(-1,3)
    for j,coefficient in enumerate(coef):
        iz,iy,ix=np.unravel_index(j,(box,)*3);center=(np.array([ix,iy,iz])+.5)/box-.5
        xyz=offset+center;phase=np.exp(2j*np.pi*frequency@xyz.T)
        for r,beta in enumerate(MONOMIALS):
            expected[:,r]+=coefficient*box**1.5*(phase@(weights*np.prod(xyz**beta,axis=1)))
    actual=cell_fourier_moments(k,coef,box)
    np.testing.assert_allclose(actual.reshape(6,10),expected,rtol=2e-9,atol=2e-11)
    np.testing.assert_allclose(direct_cell_fourier_moments(k,coef,box).reshape(6,10),expected,rtol=2e-9,atol=2e-11)
    observed=cell_forward(k,np.ones((2,3)),coef,box,1.).reshape(2,6)
    np.testing.assert_allclose(actual[...,0],observed[:,:3]-1j*observed[:,3:],atol=2e-10,rtol=1e-9)


def test_known_pilot_pose_pairing_against_nonlinear_cell_finite_differences():
    from fourier_splats.uq_continuous_pose import pose_cell_forward,PAIRS,PAIR_SCALE
    rng=np.random.default_rng(610292);box=4;n=2;nq=3
    k=rng.normal(size=(n,nq,3));q=rng.normal(size=(n,nq,2));ctf=rng.normal(size=(n,nq))
    w=rng.normal(size=2*n*nq);pilot=rng.normal(size=box**3);pilot/=np.linalg.norm(pilot)
    op=PolynomialPoseFieldOperator(k,q,ctf,w,.8,.04,.003,order=12,backend='direct')
    pair=pair_pose_fourier_moments(op,cell_fourier_moments(k,pilot,box))
    assert check_cell_pose_pairings(op,pilot,box,pair)['passed']
    u=rng.normal(size=(n,5));u/=np.linalg.norm(u,axis=1)[:,None]
    direction=np.r_[u.ravel(),np.zeros(15*n)]
    step=1e-3
    plus=w@pose_cell_forward(k,q,ctf,pilot,box,.8,step*u,.04,.003)
    minus=w@pose_cell_forward(k,q,ctf,pilot,box,.8,-step*u,.04,.003)
    zero=w@pose_cell_forward(k,q,ctf,pilot,box,.8,np.zeros_like(u),.04,.003)
    np.testing.assert_allclose(pair@direction,(plus-minus)/(2*step),atol=2e-9,rtol=2e-7)
    tensor=np.array([[v[a]*v[b]*s for (a,b),s in zip(PAIRS,PAIR_SCALE)] for v in u])
    second=np.r_[np.zeros(5*n),tensor.ravel()]
    np.testing.assert_allclose(2*pair@second,(plus+minus-2*zero)/step**2,atol=2e-8,rtol=2e-5)


def test_pilot_specific_pairing_bound_retains_full_density_class():
    rng=np.random.default_rng(610293);F=rng.normal(size=(8,4));h=rng.normal(size=8);pilot=rng.normal(size=8)
    L=1.2;B=2.;P=np.linalg.norm(pilot);f=L*np.linalg.norm(F,2)
    original=joint_density_pose_bias(np.linalg.norm(h),f,np.linalg.norm(F.T@h),L,B,P,0.)
    refined=joint_density_pose_bias(np.linalg.norm(h),f,np.linalg.norm(F.T@h),L,B,P,0.,np.linalg.norm(F.T@pilot))
    assert refined['bias_upper']<=original['bias_upper']
    for v in rng.normal(size=(100,4)):
        v*=L/np.linalg.norm(v)
        assert B*np.linalg.norm(h-F@v)+abs(pilot@(F@v))<=refined['bias_upper']+1e-10
