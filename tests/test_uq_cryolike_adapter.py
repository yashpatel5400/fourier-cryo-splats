import numpy as np
import pytest
pytest.importorskip('cryolike',reason='Optional pinned author code must be on PYTHONPATH')
import torch
from scipy.spatial.transform import Rotation
from cryolike.metadata import ViewingAngles
from fourier_splats.uq_cryolike_adapter import make_grid,prepare_images,prepare_templates


def test_author_transforms_against_independent_asymmetric_direct_sums():
    torch.set_num_threads(1)
    rng=np.random.default_rng(918734);n=12;pixel=2.3
    grid=make_grid(radius=4,inplanes=16,radial_step=.5)
    image=rng.normal(size=(1,n,n))
    actual=prepare_images(image,pixel,grid,normalize=False).images_fourier.numpy().reshape(1,-1)
    y,x=np.meshgrid(np.arange(n)-n/2,np.arange(n)-n/2,indexing='ij')
    frequencies=2*np.column_stack([grid.x_points,grid.y_points])/n
    expected=(np.exp(-2j*np.pi*(frequencies@np.array([x.ravel(),y.ravel()])))@image.ravel())*4*np.pi/n
    np.testing.assert_allclose(actual[0],expected,rtol=2e-11,atol=2e-10)
    angles=np.array([[.3,.8,.2],[1.1,.7,.9]])
    views=ViewingAngles(*(torch.tensor(angles[:,i],dtype=torch.float64) for i in range(3)))
    density=rng.normal(size=(n,n,n))
    actual=prepare_templates(density,pixel,grid,views,normalize=False)
    np.testing.assert_allclose(actual.box_size,[pixel*n]*2)
    z,y,x=np.meshgrid(np.arange(n)-n/2,np.arange(n)-n/2,np.arange(n)-n/2,indexing='ij')
    xyz=np.array([x.ravel(),y.ravel(),z.ravel()])
    plane=np.column_stack([grid.x_points,grid.y_points,np.zeros(grid.n_points)])
    offset=np.sinc(2*grid.x_points)*np.sinc(2*grid.y_points)
    for i,a in enumerate(angles):
        rotation=Rotation.from_euler('ZYZ',a).as_matrix()
        k=plane@rotation.T*2/n
        exact=np.exp(-2j*np.pi*(k@xyz))@density.ravel()
        exact=(exact-offset*density.sum())*np.pi/n
        np.testing.assert_allclose(actual.images_fourier.numpy()[i].ravel(),exact,rtol=5e-11,atol=1e-9)


def test_original_scoring_api_and_explicit_wall_limit():
    from cryolike.microscopy import CTF
    from fourier_splats.uq_cryolike_adapter import score
    torch.set_num_threads(1);rng=np.random.default_rng(337351)
    grid=make_grid(radius=4,inplanes=16,radial_step=.5)
    views=ViewingAngles(torch.tensor([.2,.8]),torch.tensor([.4,1.]),None)
    images=prepare_images(rng.normal(size=(3,12,12)),2.,grid)
    templates=prepare_templates(rng.normal(size=(12,12,12)),2.,grid,views)
    ctf=CTF(torch.ones((3,grid.n_shells,grid.n_inplanes),dtype=torch.float64))
    result=score(images,templates,ctf,shift_points=3,image_batch=2,template_batch=1)
    assert all(np.isfinite(v).all() and v.shape==(3,) for v in result.values())
    assert np.all(np.abs(result['cross_correlation_M'])<=1+1e-10)
    with pytest.raises(TimeoutError):
        score(images,templates,ctf,shift_points=1,wall_seconds=0.)
