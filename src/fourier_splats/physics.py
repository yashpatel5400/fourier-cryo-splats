"""Fourier convention: exp(-2 pi i k.x); cryoDRGN pose import conventions."""
import numpy as np
from scipy.fft import fftn, ifftn, fftshift, ifftshift

def fft_center(x):
    return fftshift(fftn(ifftshift(x,axes=(-2,-1)),axes=(-2,-1)),axes=(-2,-1))

def ctf(q,params):
    """q: (P,2) inverse angstroms; params: (N,9), cryoDRGN convention."""
    p=np.asarray(params,dtype=np.float64)
    dfu,dfv,angle,volt,cs,w,phase=[p[:,i,None] for i in range(2,9)]
    voltage=volt*1000
    lam=12.2639/np.sqrt(voltage+0.97845e-6*voltage**2)
    s2=(q*q).sum(-1)[None,:]; az=np.arctan2(q[:,1],q[:,0])[None,:]
    df=0.5*(dfu+dfv+(dfu-dfv)*np.cos(2*(az-angle*np.pi/180)))
    gamma=2*np.pi*(-0.5*df*lam*s2+0.25*cs*1e7*lam**3*s2*s2)-phase*np.pi/180
    return (np.sqrt(1-w*w)*np.sin(gamma)-w*np.cos(gamma)).astype(np.float32)

def plane_restriction(mu,precision,E):
    """Exact restriction, not a marginal/integrated projection."""
    A=E.T@precision@E; b=E.T@precision@mu
    m=np.linalg.solve(A,b)
    amplitude=np.exp(-0.5*(mu@precision@mu-b@m))
    return m,A,amplitude

def hermitian_pairs(k,centers,precisions,coefficients):
    d1=k[:,None,:]-centers[None,:,:];d2=k[:,None,:]+centers[None,:,:]
    a=np.exp(-0.5*np.einsum('nki,kij,nkj->nk',d1,precisions,d1))
    b=np.exp(-0.5*np.einsum('nki,kij,nkj->nk',d2,precisions,d2))
    return a@coefficients+b@np.conjugate(coefficients)

def volume_from_fourier(f):
    return fftshift(ifftn(ifftshift(f))).real.astype(np.float32)
