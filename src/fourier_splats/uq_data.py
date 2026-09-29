"""Group-preserving particle access and an independent voxel/NUFFT generator."""
import csv
import json
import sys
from pathlib import Path
import finufft
import mrcfile
import numpy as np
from scipy.fft import fftn, ifftn, fftshift, ifftshift
from scipy.sparse.linalg import LinearOperator
from .physics import ctf, fft_center
from .uq_physics import gaussian_pair_features


class VoxelObservationOperator:
    """Real Fourier-slice forward/adjoint pair for a declared voxel model.

    The density parameter is the vector of voxel values, with the ordinary
    Euclidean norm. Physical L2 energy additionally includes voxel volume.
    A wider voxel model can audit weights obtained from a Gaussian dictionary:
    computing ||ell-A* w|| here charges unresolved density outside that
    dictionary instead of assuming it away. It still declares a finite grid,
    fixed poses/CTFs and Gaussian whitened noise.
    """
    def __init__(self,k,transfer,box,noise_std=1.0,eps=1e-11):
        self.k=np.asarray(k,dtype=float);self.box=int(box);self.eps=eps
        if self.k.ndim!=3 or self.k.shape[-1]!=3 or self.box<2 or self.box%2:
            raise ValueError('Even grid and particle-by-frequency-by-3 coordinates required')
        self.n,self.q=self.k.shape[:2]
        self.transfer=np.broadcast_to(np.asarray(transfer,dtype=float),self.k.shape[:2])/noise_std
        if noise_std<=0 or not np.isfinite(self.transfer).all():raise ValueError('Positive finite noise required')
        xyz=2*np.pi*self.k.reshape(-1,3)/self.box
        self.coordinates=[np.ascontiguousarray(xyz[:,i]) for i in [2,1,0]]
        self.shape=(self.n*2*self.q,self.box**3)
        self.nthreads=1 if sys.platform=='darwin' else 4

    def forward(self,volume):
        v=np.asarray(volume).reshape((self.box,)*3).astype(complex)
        f=finufft.nufft3d2(*self.coordinates,v,isign=-1,eps=self.eps,nthreads=self.nthreads)
        f=f.reshape(self.n,self.q)*self.transfer
        return np.concatenate([f.real,f.imag],axis=1).ravel()

    def adjoint(self,weights):
        w=np.asarray(weights,dtype=float).reshape(self.n,2*self.q)
        c=np.ascontiguousarray(((w[:,:self.q]+1j*w[:,self.q:])*self.transfer).ravel())
        v=finufft.nufft3d1(*self.coordinates,c,(self.box,)*3,isign=1,eps=self.eps,nthreads=self.nthreads)
        return v.real.ravel()

    def forward_columns(self,columns,batch=16):
        columns=np.asarray(columns,dtype=float)
        if columns.shape[0]!=self.shape[1]:raise ValueError('Voxel-column dimensions do not match')
        result=[]
        for start in range(0,columns.shape[1],batch):
            v=np.ascontiguousarray(columns[:,start:start+batch].T.reshape((-1,self.box,self.box,self.box)),dtype=complex)
            f=finufft.nufft3d2(*self.coordinates,v,isign=-1,eps=self.eps,nthreads=self.nthreads)
            f=f.reshape(len(v),self.n,self.q)*self.transfer[None]
            result.append(np.concatenate([f.real,f.imag],axis=2).reshape(len(v),-1).T)
        return np.concatenate(result,axis=1)

    def bias_residual(self,weights,functional):
        return np.asarray(functional).ravel()-self.adjoint(weights)


class SupportedVoxelOperator:
    """A voxel operator parameterized only by the explicitly allowed support."""
    def __init__(self,operator,mask):
        self.op=operator;mask=np.asarray(mask,dtype=bool).ravel()
        if len(mask)!=operator.shape[1] or not mask.any():raise ValueError('Nonempty compatible support required')
        self.active=np.flatnonzero(mask);self.shape=(operator.shape[0],len(self.active))
        z,y,x=np.unravel_index(self.active,(operator.box,)*3)
        self.xyz=np.stack([x,y,z],axis=1)-operator.box//2
        self.gram_diagonal=np.full(len(self.active),np.sum(operator.transfer**2))

    def expand(self,values):
        out=np.zeros(self.op.shape[1]);out[self.active]=np.asarray(values).ravel();return out

    def forward(self,values):return self.op.forward(self.expand(values))
    def adjoint(self,weights):return self.op.adjoint(weights)[self.active]

    def forward_columns(self,columns,batch=16):
        pieces=[]
        for start in range(0,columns.shape[1],batch):
            block=np.zeros((self.op.shape[1],min(batch,columns.shape[1]-start)))
            block[self.active]=columns[:,start:start+batch]
            pieces.append(self.op.forward_columns(block,batch))
        return np.concatenate(pieces,axis=1)

    def as_linear_operator(self):
        return LinearOperator(self.shape,matvec=self.forward,rmatvec=self.adjoint,dtype=float)


def voxel_pose_terms(operator,pilot,q,rotation_radius,shift_radius,density_radius,order=1):
    """Ambient-space pose Jacobian, interaction Gram factor and uniform remainder.

    The density norm is Euclidean on the declared support, and poses use the
    same normalized five-dimensional rotation/shift ball as uq_physics. Full
    complex exponential columns give analytic derivatives without a Gaussian
    representation restriction. Only pixel-space interaction Grams are stored.
    CTFs and whitening remain fixed. Bounds use the full pose segment, not a
    nominal Hessian alone. Exact guarantees concern real arithmetic.
    """
    if not 0<=rotation_radius<=np.pi or shift_radius<0 or density_radius<0:
        raise ValueError('Invalid density/pose radii')
    if order not in [1,2]:raise ValueError('Supported Taylor orders are one and two')
    op=operator.op;pilot=np.asarray(pilot);q=np.asarray(q)
    if pilot.shape!=(operator.shape[1],) or q.shape!=(op.n,op.q,2):raise ValueError('Incompatible pilot/frequencies')
    xyz=operator.xyz;normx=np.linalg.norm(xyz,axis=1)
    factors=[];jacobians=[];remainder=[];pilot_bounds=[];operator_bounds=[]
    quadratic_jacobians=[];quadratic_factors=[];third_bounds=[]
    scale=-2*np.pi/op.box
    for i in range(op.n):
        phase=np.exp(1j*scale*(op.k[i]@xyz.T))*op.transfer[i,:,None]
        rotation=scale*rotation_radius*np.cross(xyz[None,:,:],op.k[i,:,None,:])
        translation=np.broadcast_to(scale*shift_radius*q[i,:,None,:],(op.q,len(xyz),2))
        phase_first=np.concatenate([rotation,translation],axis=-1)
        derivative=1j*phase[...,None]*phase_first
        d=np.concatenate([derivative.real,derivative.imag],axis=0)
        jacobians.append(np.einsum('mpq,p->mq',d,pilot))
        flat=d.reshape(2*op.q,-1);gram=flat@flat.T;maximum=np.max(np.diag(gram))
        factors.append(np.linalg.cholesky(gram+1e-12*maximum*np.eye(2*op.q)) if maximum else np.zeros_like(gram))
        kr=np.linalg.norm(op.k[i],axis=1)
        speed=abs(scale)*(rotation_radius*kr[:,None]*normx[None,:]+shift_radius*np.linalg.norm(q[i],axis=1)[:,None])
        acceleration=abs(scale)*rotation_radius**2*kr[:,None]*normx[None,:]
        second=np.abs(op.transfer[i,:,None])*(speed**2+acceleration)
        h=np.linalg.norm(second@np.abs(pilot));k=np.linalg.norm(second)
        pilot_bounds.append(h);operator_bounds.append(k)
        if order==1:
            remainder.append(.5*(h+density_radius*k))
        else:
            phase_second=np.zeros((*phase.shape,5,5))
            dot=op.k[i]@xyz.T
            for a in range(3):
                for b in range(3):
                    phase_second[...,a,b]=scale*rotation_radius**2*(
                        .5*(op.k[i,:,a,None]*xyz[None,:,b]+op.k[i,:,b,None]*xyz[None,:,a])-(a==b)*dot)
            second_derivative=phase[...,None,None]*(1j*phase_second-phase_first[..., :,None]*phase_first[...,None,:])
            d2=np.concatenate([second_derivative.real,second_derivative.imag],axis=0)
            quadratic_jacobians.append(np.einsum('mpab,p->mab',d2,pilot))
            flat2=d2.reshape(2*op.q,-1);gram2=flat2@flat2.T;maximum2=np.max(np.diag(gram2))
            quadratic_factors.append(np.linalg.cholesky(gram2+1e-12*maximum2*np.eye(2*op.q)) if maximum2 else np.zeros_like(gram2))
            jerk=abs(scale)*rotation_radius**3*kr[:,None]*normx[None,:]
            third=np.abs(op.transfer[i,:,None])*(speed**3+3*speed*acceleration+jerk)
            h3=np.linalg.norm(third@np.abs(pilot));k3=np.linalg.norm(third)
            third_bounds.append((h3,k3));remainder.append((h3+density_radius*k3)/6)
    result={'jacobian':np.asarray(jacobians),'interaction_factor':np.asarray(factors),
            'remainder':np.asarray(remainder),'pilot_second_bound':np.asarray(pilot_bounds),
            'operator_second_bound':np.asarray(operator_bounds),'taylor_order':order}
    if order==2:
        result.update(quadratic_jacobian=np.asarray(quadratic_jacobians),quadratic_interaction_factor=np.asarray(quadratic_factors),
                      pilot_third_bound=np.asarray(third_bounds)[:,0],operator_third_bound=np.asarray(third_bounds)[:,1])
    return result


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
