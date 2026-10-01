import numpy as np
from fourier_splats.uq_shift_moments import shifted_real_means,gaussian_shift_second_moment


def test_zero_shift_matches_dense_real_second_moment():
    rng=np.random.default_rng(151);f=rng.normal(size=(7,4))+1j*rng.normal(size=(7,4))
    q=np.array([[1,0],[0,1],[2,3],[4,-1]])
    real=np.concatenate([f.real,f.imag],axis=1)
    np.testing.assert_allclose(gaussian_shift_second_moment(f,q,0),real.T@real/len(real),atol=1e-14)


def test_gaussian_moment_matches_independent_quadrature():
    from numpy.polynomial.hermite import hermgauss
    f=np.array([[1+2j,-.3+1j,2-.1j],[-1+.3j,.8+.2j,1.4-1j]])
    q=np.array([[1,0],[0,2],[3,-1]]);nodes,weights=hermgauss(32)
    xx,yy=np.meshgrid(nodes,nodes);shifts=np.sqrt(2)*2*np.stack([xx.ravel(),yy.ravel()],axis=1)
    joint=(weights[:,None]*weights[None,:]).ravel()/np.pi;moment=np.zeros((6,6))
    for frame in f:
        means=shifted_real_means(np.tile(frame,(len(shifts),1)),q,shifts,box=16)
        moment+=np.einsum('n,ni,nj->ij',joint,means,means)/len(f)
    np.testing.assert_allclose(gaussian_shift_second_moment(f,q,2,box=16),moment,rtol=1e-11,atol=1e-11)


def test_broad_shift_removes_cross_frequency_and_pseudo_moments():
    f=np.array([[1+2j,3-4j]]);q=np.array([[1,0],[0,2]])
    actual=gaussian_shift_second_moment(f,q,1e4)
    np.testing.assert_allclose(actual,np.diag([2.5,12.5,2.5,12.5]),atol=1e-12)
