"""Sparse, conjugate-tied spectral Gaussian and trilinear operators.

The Gaussian model uses fixed lattice centers and fixed isotropic width. It is
an explicit stationary instance of the general paired Gaussian representation.
"""
import itertools
import numpy as np
from scipy.sparse import csr_matrix

def design(k,box,kind='gaussian',sigma=0.5,radius=2.0):
    k=np.asarray(k,dtype=np.float32)
    half=box//2; side=box+1; mid=(side**3-1)//2; ncoef=mid+1
    if kind=='gaussian':
        extent=int(np.ceil(radius))
        offsets=np.array(list(itertools.product(range(-extent,extent+1),repeat=3)),dtype=np.int32)
        base=np.rint(k).astype(np.int32)
    elif kind=='voxel':
        offsets=np.array(list(itertools.product([0,1],repeat=3)),dtype=np.int32)
        base=np.floor(k).astype(np.int32)
    else: raise ValueError(kind)
    # xyz coordinates, linear indices stored in zyx volume order.
    nodes=base[:,None,:]+offsets[None,:,:]
    delta=k[:,None,:]-nodes
    if kind=='gaussian':
        d2=np.sum(delta**2,axis=-1)
        weights=np.exp(-0.5*d2/sigma**2).astype(np.float32)
        keep=d2<=radius**2
    else:
        weights=np.prod(np.maximum(1-np.abs(delta),0),axis=-1).astype(np.float32)
        keep=weights>0
    keep &= (np.abs(nodes)<=half).all(-1)
    row,col=np.nonzero(keep)
    nodes=nodes[row,col]+half; val=weights[row,col]
    linear=(nodes[:,2]*side+nodes[:,1])*side+nodes[:,0]
    inverse=side**3-1-linear
    canonical=np.minimum(linear,inverse)
    sign=np.where(linear<=mid,1.,-1.).astype(np.float32)
    sign[linear==mid]=0 # DC must be real.
    ar=csr_matrix((val,(row,canonical)),shape=(len(k),ncoef),dtype=np.float32)
    ai=csr_matrix((val*sign,(row,canonical)),shape=(len(k),ncoef),dtype=np.float32)
    ar.sum_duplicates();ai.sum_duplicates();ai.eliminate_zeros()
    return ar,ai

def evaluate_grid(real,imag,box,kind='gaussian',sigma=0.5,radius=2.0):
    q=np.arange(-box//2,box//2,dtype=np.float32)
    z,y,x=np.meshgrid(q,q,q,indexing='ij');k=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    f=np.empty(len(k),np.complex64)
    for start in range(0,len(k),16384):
        ar,ai=design(k[start:start+16384],box,kind,sigma,radius)
        f[start:start+len(ar.indptr)-1]=ar@real+1j*(ai@imag)
    f=f.reshape(box,box,box)
    # The even FFT grid has self-paired Nyquist planes. Exclude these and the
    # unsupported outer sphere; do not alter any interior Fourier phase.
    f[0]=0;f[:,0]=0;f[:,:,0]=0
    f[(x*x+y*y+z*z)>(box/2-2)**2]=0
    return f
