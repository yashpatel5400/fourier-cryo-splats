"""Analytic Fourier moments of a known orthonormal constant-cell pilot.

The pilot is known; this does not restrict the unknown continuous density.
NUFFT/special-function arithmetic is numerically checked, not validated by
interval arithmetic. The integrals themselves include each cell's extent.
"""
import itertools
import numpy as np
import finufft
from .uq_joint_bias import MONOMIALS


def normalized_cell_moments(k,box):
    """E[X^j exp(2 pi i k X)] for X uniform on [-1/(2m),1/(2m)]."""
    u=np.asarray(k,float)/box;x=np.pi*u
    first=np.empty_like(x);second=np.empty_like(x);small=abs(x)<.03
    t=x[small];t2=t*t
    first[small]=np.pi*t*(-1/3+t2*(1/30+t2*(-1/840+t2*(1/45360-t2/3991680))))
    second[small]=np.pi**2*(-1/3+t2*(1/10+t2*(-1/168+t2*(1/6480-t2/443520))))
    t=x[~small]
    first[~small]=np.pi*(t*np.cos(t)-np.sin(t))/t**2
    second[~small]=np.pi**2*(-np.sin(t)/t-2*np.cos(t)/t**2+2*np.sin(t)/t**3)
    return np.stack([np.sinc(u),first/(2j*np.pi*box),-second/(4*np.pi**2*box**2)],axis=-1)


def cell_fourier_moments(k,coefficients,box,eps=1e-12,nthreads=1):
    """Integral rho0(x) x**beta exp(+2 pi i k.x) for ten degree<=2 moments.

    Cell coefficients use the orthonormal basis from cell_forward: density is
    m**1.5 times the coefficient on each cell, with array order (z,y,x).
    Expand x=center+offset, transform center monomials and multiply analytic
    offset moments. All Fourier signs and the half-cell phase are explicit.
    """
    k=np.asarray(k,float);coefficients=np.asarray(coefficients,float)
    if box<2 or box%2 or coefficients.size!=box**3 or k.shape[-1]!=3:
        raise ValueError('Even cell grid, matching coefficients and 3D frequencies required')
    if not np.isfinite(k).all() or not np.isfinite(coefficients).all():raise ValueError('Finite inputs required')
    coordinate=(np.arange(box)+.5)/box-.5
    z,y,x=np.meshgrid(coordinate,coordinate,coordinate,indexing='ij');xyz=[x,y,z]
    strengths=np.stack([coefficients.reshape((box,)*3)*np.prod(np.stack([xyz[a]**power for a,power in enumerate(beta)]),axis=0) for beta in MONOMIALS]).astype(complex)
    frequencies=[np.ascontiguousarray(2*np.pi*k.reshape(-1,3)[:,axis]/box) for axis in [2,1,0]]
    transformed=finufft.nufft3d2(*frequencies,strengths,isign=1,eps=eps,nthreads=nthreads)
    transformed=transformed.T.reshape((*k.shape[:-1],10))*np.exp(1j*np.pi*k.sum(axis=-1)/box)[...,None]/box**1.5
    local=normalized_cell_moments(k,box);lookup={powers:j for j,powers in enumerate(MONOMIALS)}
    result=np.zeros_like(transformed)
    for j,beta in enumerate(MONOMIALS):
        for gamma in itertools.product(*(range(p+1) for p in beta)):
            multiplier=np.ones(k.shape[:-1],complex)
            for axis,(b,g) in enumerate(zip(beta,gamma)):
                multiplier*=local[...,axis,b-g]*(2 if b==2 and g==1 else 1)
            result[...,j]+=transformed[...,lookup[gamma]]*multiplier
    return result


def direct_cell_fourier_moments(k,coefficients,box,chunk=32):
    """Direct sums over physical cell centers, without a NUFFT/grid phase.

    This numerical cross-check deliberately forms each cell's translated
    moments before summing. It shares only the scalar within-cell formulas
    (which have their own independent quadrature tests) with the fast path.
    """
    k=np.asarray(k,float);coef=np.asarray(coefficients,float).reshape(-1)
    if coef.size!=box**3:raise ValueError('Matching constant-cell coefficients required')
    centers=(np.indices((box,)*3).reshape(3,-1).T[:,::-1]+.5)/box-.5
    frequencies=k.reshape(-1,3);result=np.empty((len(frequencies),10),complex)
    for start in range(0,len(frequencies),chunk):
        frequency=frequencies[start:start+chunk]
        phase=np.exp(2j*np.pi*frequency@centers.T)
        local=normalized_cell_moments(frequency,box)
        translated=[]
        for axis in range(3):
            c=centers[None,:,axis];v=local[:,axis,:]
            translated.append((v[:,0,None],c*v[:,0,None]+v[:,1,None],
                               c*c*v[:,0,None]+2*c*v[:,1,None]+v[:,2,None]))
        for j,beta in enumerate(MONOMIALS):
            moment=phase.copy()
            for axis,power in enumerate(beta):moment*=translated[axis][power]
            result[start:start+len(frequency),j]=moment@coef/box**1.5
    return result.reshape((*k.shape[:-1],10))


def check_cell_pose_pairings(op,coefficients,box,pairing,moments=None):
    """Check every lifted column for three deterministic archive particles.

    This is an independent exponential-sum diagnostic, NOT an all-column
    numerical upper bound or an interval-arithmetic error certificate.
    """
    from copy import copy
    from .uq_joint_bias import pair_pose_fourier_moments
    indices=np.unique([0,op.n//2,op.n-1]);subset=copy(op);subset.n=len(indices)
    for key in ['k','q','c','angle','shift']:setattr(subset,key,getattr(op,key)[indices])
    columns=np.r_[(5*indices[:,None]+np.arange(5)).ravel(),
                  (5*op.n+15*indices[:,None]+np.arange(15)).ravel()]
    subset.denominators=op.denominators[columns]
    direct=direct_cell_fourier_moments(subset.k,coefficients,box)
    if moments is None:moments=cell_fourier_moments(subset.k,coefficients,box,nthreads=op.nthreads)
    else:moments=np.asarray(moments)[indices]
    expected=pair_pose_fourier_moments(subset,direct);actual=np.asarray(pairing)[columns]
    absolute=float(np.max(abs(expected-actual)))
    relative=float(np.linalg.norm(expected-actual)/max(np.linalg.norm(expected),1e-300))
    moments_relative=float(np.linalg.norm(moments-direct)/max(np.linalg.norm(direct),1e-300))
    return {'particle_positions':indices.tolist(),'lifted_column_count':len(columns),
            'maximum_absolute_pairing_difference':absolute,'relative_pairing_difference':relative,
            'relative_moment_difference':moments_relative,
            'passed':bool(absolute<1e-12 or relative<1e-8),
            'scope':'Selected columns only; direct physical-cell exponential sums check NUFFT, axis order and half-cell phase. Scalar within-cell formulas shared; not a global floating-point pad.'}
