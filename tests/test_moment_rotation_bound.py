import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.spatial.transform import Rotation
from fourier_splats.uq_moment_rotation_bound import (
    cell_absolute_radial_moments,contrast_rotation_curvature,euler_cover_for_remainder)
from fourier_splats.uq_bispectrum import moment_features


def test_physical_cell_radial_integrals_include_cell_extent():
    # A single full-cube cell has nonzero radial moments at center zero.
    bound=cell_absolute_radial_moments(np.ones(1),1)
    np.testing.assert_allclose(bound,[1.,.5,.25])
    x,w=leggauss(24);x=x/2;w=w/2
    z,y,x=np.meshgrid(x,x,x,indexing='ij');weights=np.einsum('i,j,k->ijk',w,w,w)
    first=np.sum(weights*np.sqrt(x*x+y*y+z*z))
    second=np.sum(weights*(x*x+y*y+z*z))
    assert 0<first<=bound[1]
    np.testing.assert_allclose(second,bound[2],atol=1e-14)


def test_geodesic_derivative_envelope_with_signed_atomic_density():
    rng=np.random.default_rng(260101)
    x=rng.uniform(-.5,.5,(17,3));density=rng.normal(size=17)/17
    radii=np.linalg.norm(x,axis=1)
    moments=np.array([np.abs(density)@radii**k for k in range(3)])
    q=np.array([[1.,0,0],[0,1,0],[1,1,0],[2,1,0]])
    triads=np.array([[0,1,2],[0,2,3]])
    transfer=np.array([.8,-1.1,.3,.7]);w=rng.normal(size=8)
    bound=contrast_rotation_curvature(q,transfer,triads,w,moments,1.1)['total_curvature']
    for rotation in Rotation.random(12,random_state=rng).as_matrix():
        axis=rng.normal(size=3);axis/=np.linalg.norm(axis)
        k=q@rotation;kp=np.cross(k,axis);kpp=np.cross(kp,axis)
        e=np.exp(-2j*np.pi*k@x.T)
        phase1=-2j*np.pi*kp@x.T;phase2=-2j*np.pi*kpp@x.T
        m=1.1*transfer*(e@density)
        first=1.1*transfer*((e*phase1)@density)
        second=1.1*transfer*((e*(phase1**2+phase2))@density)
        value=np.sum(w[:4]*(np.abs(first)**2+(second*m.conj()).real))
        for j,(a,b,c) in enumerate(triads):
            derivative=(second[a]*m[b]*m[c].conj()+m[a]*second[b]*m[c].conj()
                +m[a]*m[b]*second[c].conj()+2*(first[a]*first[b]*m[c].conj()
                +first[a]*m[b]*first[c].conj()+m[a]*first[b]*first[c].conj()))
            value+=np.real((w[4+j]-1j*w[6+j])*derivative)/2
        assert abs(value)<=bound
        eps=1e-4
        def f(t):
            kk=q@rotation@Rotation.from_rotvec(t*axis).as_matrix()
            z=1.1*transfer*(np.exp(-2j*np.pi*kk@x.T)@density)
            return (moment_features(z,triads)@w).item()
        np.testing.assert_allclose((f(eps)-2*f(0)+f(-eps))/eps**2,value,rtol=2e-5,atol=2e-5)


def test_stationary_maximum_cover_for_known_rotation_objective():
    # f(R)=tr(A R) has exact Procrustes maximum and curvature <= nuclear norm(A).
    a=np.array([[.7,.2,-.1],[.3,.5,.4],[-.2,.1,.9]])
    u,s,vt=np.linalg.svd(a)
    exact=s[0]+s[1]+np.linalg.det(u@vt)*s[2]
    h=s.sum();grid=euler_cover_for_remainder(h,1.)
    alpha=2*np.pi*np.arange(grid['alpha_count'])/grid['alpha_count']
    beta=np.linspace(0,np.pi,grid['beta_count'])
    gamma=2*np.pi*np.arange(grid['gamma_count'])/grid['gamma_count']
    angles=np.stack(np.meshgrid(alpha,beta,gamma,indexing='ij'),axis=-1).reshape(-1,3)
    rotations=Rotation.from_euler('ZYZ',angles).as_matrix()
    sampled=np.einsum('ij,nji->n',a,rotations).max()
    assert sampled<=exact+1e-12
    assert exact<=sampled+grid['remainder_bound']+1e-12
    assert grid['points']==len(rotations) and grid['remainder_bound']<=1.+1e-12


def test_cover_count_and_tolerance():
    for h in [0.,.01,5.,10000.]:
        for epsilon in [.1,.001]:
            grid=euler_cover_for_remainder(h,epsilon)
            assert grid['remainder_bound']<=epsilon*(1+1e-12)
            assert grid['points']>=1
