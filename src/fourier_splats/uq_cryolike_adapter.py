"""Thin CPU adapter to pinned author CryoLike; no replacement scoring kernel.

Project images are [image,y,x], volumes [z,y,x]. CryoLike's physical arrays
are [image,x,y] and [x,y,z]. Its polar radius is half the integer-box Fourier
radius. Explicit conversions are tested against direct Fourier sums.
"""
import numpy as np
import torch
import time
from cryolike.grids import PolarGrid, PhysicalImages, PhysicalVolume, Volume
from cryolike.stacks import Images, Templates
from cryolike.util import Precision
from cryolike.microscopy import CTF
from cryolike.likelihoods import template_first_comparator, compute_optimal_pose
from .physics import ctf


def make_grid(radius=12., inplanes=128, radial_step=.25):
    return PolarGrid(radius_max=radius/2,dist_radii=radial_step,
                     n_inplanes=inplanes,uniform=True)


def prepare_images(images_yx,pixel_size,grid,normalize=True):
    values=np.array(images_yx,dtype=np.float64,copy=True).transpose(0,2,1).copy()
    images=Images(phys_data=PhysicalImages(torch.from_numpy(values),pixel_size))
    if normalize:images.center_physical_image_signal()
    images.transform_to_fourier(grid,precision=Precision.DOUBLE,device='cpu',nufft_eps=1e-12)
    if normalize:images.normalize_images_fourier(ord=2,use_max=False)
    return images


def prepare_templates(volume_zyx,pixel_size,grid,views,normalize=True):
    values=np.array(volume_zyx,dtype=np.float64,copy=True).transpose(2,1,0).copy()
    volume=Volume(density_physical_data=PhysicalVolume(torch.from_numpy(values),pixel_size))
    templates=Templates.generate_from_physical_volume(volume,grid,views,
        precision=Precision.DOUBLE,device=torch.device('cpu'),nufft_eps=1e-12)
    # The low-level constructor returns the normalized default 2x2 box;
    # displacement kernels need the actual field width in angstroms.
    templates.box_size=np.array([values.shape[0]*pixel_size]*2)
    if normalize:templates.normalize_images_fourier(ord=2,use_max=False)
    return templates


def prepare_ctf(parameters,field_A,grid):
    q=2*np.column_stack([grid.x_points,grid.y_points])/field_A
    values=ctf(q,parameters).astype(np.float64)
    return CTF(torch.from_numpy(values.reshape(len(parameters),grid.n_shells,grid.n_inplanes)))


def score(images,templates,transfer,*,max_shift_pixels=2.,shift_points=5,
          image_batch=16,template_batch=8,wall_seconds=900.,callback=None):
    templates.set_displacement_grid(max_shift_pixels,shift_points,shift_points,
                                   pixel_size_angstrom=float(images.phys_grid.pixel_size[0]))
    iterator=template_first_comparator(torch.device('cpu'),images,templates,transfer,
        n_images_per_batch=image_batch,n_templates_per_batch=template_batch,
        return_integrated_likelihood=True,precision=Precision.DOUBLE)
    start=time.perf_counter()
    def bounded_iterator():
        for batch in iterator:
            elapsed=time.perf_counter()-start
            if callback is not None:callback(dict(template_end=batch.t_end,image_end=batch.i_end,seconds=elapsed))
            if elapsed>wall_seconds:raise TimeoutError('Declared CryoLike scoring wall limit')
            yield batch
    pose,integrated=compute_optimal_pose(bounded_iterator(),templates,images,Precision.DOUBLE,True)
    return {**{name:val.detach().cpu().numpy() for name,val in zip(pose._fields,pose)},
            'integrated_log_score':integrated.detach().cpu().numpy()}
