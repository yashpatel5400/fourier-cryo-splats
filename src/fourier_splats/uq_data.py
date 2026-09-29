"""Group-preserving particle access and an independent voxel/NUFFT generator."""
import csv
import json
import sys
from pathlib import Path
import finufft
import mrcfile
import numpy as np
from scipy.fft import fftn, ifftn, fftshift, ifftshift
from .physics import ctf, fft_center
from .uq_physics import gaussian_pair_features


def half_plane_frequencies(radius):
    limit=int(np.ceil(radius))
    return np.array([(x,y) for y in range(-limit,limit+1) for x in range(-limit,limit+1)
                     if 0<x*x+y*y<=radius**2 and (y>0 or (y==0 and x>0))],dtype=float)


def particle_geometry(root,dataset,split,radius=5,count=None,seed=0):
    root=Path(root);directory=root/'data'/str(dataset)
    rows=list(csv.DictReader((root/'research/uncertainty/splits'/(str(dataset)+'.csv')).open()))
    valid=[r for r in rows if r['split'] in ([split] if isinstance(split,str) else split)]
    indices=np.array([int(r['output_index']) for r in valid])
    groups=np.array([r['source_group'] for r in valid])
    if count is not None and count<len(indices):
        chosen=np.random.default_rng(seed).choice(len(indices),count,replace=False)
        indices,groups=indices[chosen],groups[chosen]
    metadata=np.load(directory/'metadata.npz')
    manifest=json.loads((directory/'manifest.json').read_text())
    q=half_plane_frequencies(radius)
    plane=np.pad(q,((0,0),(0,1)))
    rotations=metadata['rotations'][indices]
    k=np.einsum('qi,nij->nqj',plane,rotations)
    field=manifest['raw_box']*manifest['raw_pixel_size_A']
    transfer=ctf(q/field,metadata['ctf'][indices]).astype(float)
    return {'indices':indices,'groups':groups,'q':np.broadcast_to(q,(len(indices),len(q),2)).copy(),
            'k':k,'ctf':transfer,'field_A':field,'manifest':manifest,
            'translations':metadata['translations'][indices]}


def particle_observations(root,dataset,geometry):
    """Centered image Fourier samples, without data-dependent normalization.

    The supplied published shifts center the observations; inference is
    conditional on these externally fitted shifts and rotations.
    """
    images=np.load(Path(root)/'data'/str(dataset)/'images.npy',mmap_mode='r')
    selected=np.asarray(images[geometry['indices']],dtype=float)
    fourier=fft_center(selected)
    q=geometry['q'][0].astype(int);center=selected.shape[-1]//2
    samples=fourier[:,q[:,1]+center,q[:,0]+center]*geometry['manifest']['data_sign']
    samples*=np.exp(-2j*np.pi*np.einsum('nqi,ni->nq',geometry['q'],geometry['translations']))
    return samples


def resample_volume(volume,box):
    """Low-pass/crop a cubic volume, retaining its physical field of view.

    Remove the retained even-grid Nyquist planes so the cropped spectrum has
    exact Hermitian partners. The resulting array defines the known generator;
    it is not represented as the original native-resolution map.
    """
    source=volume.shape[0]
    if volume.shape!=(source,)*3 or box>source or box%2 or source%2:
        raise ValueError('Expected even cubic source and an even smaller box')
    f=fftshift(fftn(ifftshift(np.asarray(volume,dtype=np.float32)),workers=4))
    start=source//2-box//2
    cropped=f[start:start+box,start:start+box,start:start+box].copy()*(box/source)**3
    cropped[0,:,:]=0;cropped[:,0,:]=0;cropped[:,:,0]=0
    result=fftshift(ifftn(ifftshift(cropped),workers=4))
    if np.max(np.abs(result.imag))>1e-5*max(np.max(np.abs(result.real)),1e-12):
        raise ValueError('Resampled map lost Hermitian symmetry')
    return result.real.astype(float)


class VoxelReference:
    """Known discrete voxel signal, projected without Gaussian interpolation."""
    def __init__(self,volume,field_A):
        self.volume=np.asarray(volume,dtype=float)
        self.field_A=float(field_A)
        self.box=self.volume.shape[0]
        if self.volume.shape!=(self.box,)*3:raise ValueError('Cubic volume required')

    @classmethod
    def from_mrc(cls,path,box=64):
        with mrcfile.open(path) as m:
            if tuple(int(m.header[x]) for x in ['mapc','mapr','maps'])!=(1,2,3):
                raise ValueError('Nonstandard axis mapping must be handled explicitly')
            field=float(m.voxel_size.x)*m.data.shape[2]
            volume=resample_volume(m.data,box)
        return cls(volume,field)

    def fourier(self,k,image_field_A=None):
        shape=np.asarray(k).shape[:-1]
        k=np.asarray(k,dtype=float).reshape(-1,3)
        field=self.field_A if image_field_A is None else image_field_A
        x=np.ascontiguousarray(2*np.pi*k*self.field_A/(field*self.box))
        # FINUFFT array axes are z,y,x, while supplied vectors are x,y,z.
        result=finufft.nufft3d2(np.ascontiguousarray(x[:,2]),np.ascontiguousarray(x[:,1]),
                               np.ascontiguousarray(x[:,0]),self.volume.astype(complex),
                               isign=-1,eps=1e-10,nthreads=1 if sys.platform=='darwin' else 4)
        return result.reshape(shape)


def local_weights(box,center_xyz,width_pixels):
    grid=np.arange(-box//2,box//2)
    z,y,x=np.meshgrid(grid,grid,grid,indexing='ij')
    distance=(x-center_xyz[0])**2+(y-center_xyz[1])**2+(z-center_xyz[2])**2
    weight=np.exp(-distance/(2*width_pixels**2))
    return weight/weight.sum()


def discrete_density_functionals(weights,centers,sigma):
    """Linear functionals of the real inverse-DFT Gaussian density on a grid."""
    weights=np.asarray(weights)
    if weights.ndim==3:weights=weights[None]
    box=weights.shape[-1]
    grid=np.arange(-box//2,box//2)
    z,y,x=np.meshgrid(grid,grid,grid,indexing='ij')
    k=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    f=fftshift(fftn(ifftshift(weights,axes=(-3,-2,-1)),axes=(-3,-2,-1)),axes=(-3,-2,-1))
    f=f.reshape(len(weights),-1).conjugate()/box**3
    answer=np.zeros((len(weights),2*len(centers)))
    for start in range(0,len(k),2048):
        features=gaussian_pair_features(k[start:start+2048],centers,sigma)
        answer+=(f[:,start:start+2048]@features).real
    return answer
