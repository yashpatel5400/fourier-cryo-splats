"""Exact Gaussian common-shift second moments for finite Fourier coordinates."""
import numpy as np


def shifted_real_means(fourier, q, shifts_pixels, box=64):
    f=np.asarray(fourier,complex);q=np.asarray(q,float);s=np.asarray(shifts_pixels,float)
    if f.ndim!=2 or q.shape!=(f.shape[1],2) or s.shape!=(f.shape[0],2) or box<=0:
        raise ValueError('Matching Fourier means, plane frequencies and particle shifts required')
    phase=np.exp(-2j*np.pi*(s@q.T)/box);values=f*phase
    return np.concatenate([values.real,values.imag],axis=1)


def gaussian_shift_second_moment(fourier,q,sigma_pixels,box=64):
    """Average over a uniform supplied orientation law and N(0,sigma² I) shift.

    Both conjugate and unconjugated complex moments are required for correct
    realification; keeping only the conjugate covariance loses phase content.
    The shift is common to the two exposures and independent of orientation.
    """
    f=np.asarray(fourier,complex);q=np.asarray(q,float);sigma=float(sigma_pixels)
    if f.ndim!=2 or q.shape!=(f.shape[1],2) or not np.isfinite(sigma) or sigma<0 or box<=0:
        raise ValueError('Finite matching Fourier means and nonnegative shift SD required')
    conjugate=f.T@f.conj()/len(f);pseudo=f.T@f/len(f)
    difference=np.sum((q[:,None,:]-q[None,:,:])**2,axis=-1)
    addition=np.sum((q[:,None,:]+q[None,:,:])**2,axis=-1)
    conjugate*=np.exp(-2*np.pi**2*(sigma/box)**2*difference)
    pseudo*=np.exp(-2*np.pi**2*(sigma/box)**2*addition)
    result=.5*np.block([[(conjugate+pseudo).real,(pseudo-conjugate).imag],
                       [(pseudo+conjugate).imag,(conjugate-pseudo).real]])
    return (result+result.T)/2
