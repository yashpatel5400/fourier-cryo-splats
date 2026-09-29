import numpy as np
from scipy.fft import fftn, ifftn, fftshift, ifftshift
from fourier_splats.uq_data import VoxelReference,discrete_density_functionals,local_weights
from fourier_splats.uq_physics import gaussian_pair_features


def test_nufft_generator_matches_direct_sum_and_cartesian_fft():
    rng=np.random.default_rng(204);box=8
    v=rng.normal(size=(box,)*3);ref=VoxelReference(v,field_A=40)
    q=np.array([[0.,0.,0.],[1.,2.,-1.],[-2.,1.,0.],[.21,-1.43,.78]])
    actual=ref.fourier(q)
    g=np.arange(-box//2,box//2);z,y,x=np.meshgrid(g,g,g,indexing='ij')
    xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    expected=np.exp(-2j*np.pi*q@xyz.T/box)@v.ravel()
    assert np.allclose(actual,expected,rtol=1e-8,atol=1e-8)
    f=fftshift(fftn(ifftshift(v)))
    assert np.allclose(actual[:3],f[q[:3,2].astype(int)+4,q[:3,1].astype(int)+4,q[:3,0].astype(int)+4])


def test_functionals_match_rendered_gaussian_map():
    rng=np.random.default_rng(205);box=12
    centers=rng.normal(size=(6,3));c=rng.normal(size=12)
    weights=np.stack([local_weights(box,[1,2,-1],1.4),local_weights(box,[-2,-1,2],2)])
    ell=discrete_density_functionals(weights,centers,.8)
    g=np.arange(-box//2,box//2);z,y,x=np.meshgrid(g,g,g,indexing='ij')
    k=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    f=(gaussian_pair_features(k,centers,.8)@c).reshape((box,)*3)
    rendered=fftshift(ifftn(ifftshift(f))).real
    expected=np.sum(weights*rendered,axis=(1,2,3))
    assert np.allclose(ell@c,expected,atol=1e-12)
