"""Interpolation guide for finding pose violations; never a certificate."""
import numpy as np
from scipy.ndimage import spline_filter,map_coordinates
from scipy.spatial.transform import Rotation


class InterpolatedCellFourier:
    def __init__(self,coefficients,box,nfft=256):
        coefficients=np.asarray(coefficients,float).reshape((box,)*3)
        if box%2 or nfft%2 or nfft<box:raise ValueError('Even padded grid at least as large as density required')
        padded=np.zeros((nfft,)*3);offset=(nfft-box)//2
        padded[offset:offset+box,offset:offset+box,offset:offset+box]=coefficients
        transformed=np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(padded)))/box**1.5
        self.real=spline_filter(transformed.real,order=3,mode='grid-wrap')
        self.imag=spline_filter(transformed.imag,order=3,mode='grid-wrap')
        self.box=box;self.nfft=nfft

    def values(self,k):
        k=np.asarray(k,float)
        if k.ndim!=2 or k.shape[1]!=3 or np.max(abs(k))>=self.box/2-1:
            raise ValueError('Guide queries must stay inside the native Fourier grid')
        coordinates=(k[:,::-1]*self.nfft/self.box+self.nfft/2).T
        value=map_coordinates(self.real,coordinates,order=3,mode='grid-wrap',prefilter=False)
        value=value+1j*map_coordinates(self.imag,coordinates,order=3,mode='grid-wrap',prefilter=False)
        return value*np.prod(np.sinc(k/self.box),axis=1)*np.exp(-1j*np.pi*k.sum(axis=1)/self.box)


def normalized_quadratic(fourier,q,transfer,direction,shift_pixels,box=64):
    value=np.asarray(fourier)*transfer*np.exp(-2j*np.pi*(q@shift_pixels)/box)
    mean=np.concatenate([value.real,value.imag]);denominator=mean@mean
    if denominator<=0:raise ValueError('Nonzero mean required')
    return float(mean@direction@mean/denominator)


def guide_pose_score(parameters,base_rotation,guide,q,transfer,direction):
    rotation=base_rotation@Rotation.from_rotvec(parameters[:3]).as_matrix()
    k=np.pad(q,((0,0),(0,1)))@rotation
    return normalized_quadratic(guide.values(k),q,transfer,direction,np.asarray(parameters[3:]),guide.box)
