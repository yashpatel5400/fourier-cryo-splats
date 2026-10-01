"""Analytic rotation-curvature envelopes for a fixed physical-cell density.

These are ordinary floating-point evaluations of analytic inequalities, not
outward-rounded certificates. They do not estimate experimental map error.
"""
import math
import numpy as np


def cell_absolute_radial_moments(coefficients, box):
    """Bound integral |rho| |x|^p, p=0,1,2, in the physical unit cube.

    Cell coefficients use the orthonormal box**1.5 indicator basis. M0 and M2
    are exact integrals; M1 uses Jensen inside each cell, not a center sample.
    """
    coefficients=np.asarray(coefficients,float).reshape((box,)*3)
    if box<1 or not np.all(np.isfinite(coefficients)):
        raise ValueError('Finite density and positive box required')
    x=(np.arange(box)+.5)/box-.5
    r2=x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2+1/(4*box**2)
    mass=np.abs(coefficients)/box**1.5
    return np.array([mass.sum(),np.sum(mass*np.sqrt(r2)),np.sum(mass*r2)])


def contrast_rotation_curvature(q,transfer,triads,direction,radial_moments,max_amplitude):
    """Uniform |d²/dt² f(exp(t A)R,a)| for ||axis(A)||=1, |a|<=max_amplitude.

    Triangle inequalities bound derivatives of the Fourier integral directly;
    physical cell shape is already included in the radial integrals.
    """
    q=np.asarray(q,float);transfer=np.abs(np.asarray(transfer))
    triads=np.asarray(triads,int).reshape(-1,3);direction=np.asarray(direction,float)
    m0,m1,m2=np.asarray(radial_moments,float)
    n=len(q);t=len(triads)
    if (q.ndim!=2 or q.shape[1] not in (2,3) or transfer.shape!=(n,)
        or direction.shape!=(n+2*t,) or min(m0,m1,m2,max_amplitude)<0
        or not np.all(np.isfinite(direction))):
        raise ValueError('Incompatible finite contrast or invalid nonnegative bounds')
    omega=2*np.pi*np.linalg.norm(q,axis=1)
    d0=transfer*m0;d1=transfer*omega*m1
    d2=transfer*(omega**2*m2+omega*m1)
    power=float(2*np.abs(direction[:n]/2)@(d1*d1+d0*d2))
    cubic=0.
    if t:
        a,b,c=triads.T
        weight=np.abs((direction[n:n+t]-1j*direction[n+t:])/2)
        terms=(d2[a]*d0[b]*d0[c]+d0[a]*d2[b]*d0[c]+d0[a]*d0[b]*d2[c]
            +2*(d1[a]*d1[b]*d0[c]+d1[a]*d0[b]*d1[c]+d0[a]*d1[b]*d1[c]))
        cubic=float(weight@terms)
    return dict(power_curvature=power,cubic_curvature=cubic,
        total_curvature=float(max_amplitude**2*power+max_amplitude**3*cubic))


def euler_cover_for_remainder(curvature,remainder):
    """Construct a ZYZ grid with H radius²/2 <= requested remainder.

    Count is a sufficient construction, never a lower bound on optimal cost.
    No grid evaluations are performed. Integer dimensions can be very large.
    """
    if not np.isfinite(curvature) or curvature<0 or not np.isfinite(remainder) or remainder<=0:
        raise ValueError('Finite nonnegative curvature and positive remainder required')
    if curvature==0:
        return dict(alpha_count=1,beta_count=1,gamma_count=1,points=1,
                    radius_bound=math.pi,remainder_bound=0.)
    radius=math.sqrt(2*remainder/curvature)
    n=max(1,math.ceil(3*math.pi/radius))
    m=max(2,math.ceil(3*math.pi/(2*radius))+1)
    r=min(math.pi,2*math.pi/n+math.pi/(2*(m-1)))
    return dict(alpha_count=n,beta_count=m,gamma_count=n,points=n*n*m,
                radius_bound=r,remainder_bound=curvature*r*r/2)
