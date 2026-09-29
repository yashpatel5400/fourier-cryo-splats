"""Smooth full Gaussian Fourier features for uncertainty development experiments.

Unlike the original stationary sparse implementation, these reference features
have no hard cutoff. Analytic derivatives and uniform Taylor bounds therefore
apply even when a perturbed frequency crosses a computational neighborhood.
"""
import numpy as np
from scipy.spatial.transform import Rotation


def gaussian_pair_features(k, centers, sigma):
    k, centers = np.asarray(k), np.asarray(centers)
    sigma = np.broadcast_to(np.asarray(sigma), (len(centers),))
    dplus = k[..., None, :] - centers
    dminus = k[..., None, :] + centers
    gp = np.exp(-0.5*np.sum(dplus*dplus, axis=-1)/sigma**2)
    gm = np.exp(-0.5*np.sum(dminus*dminus, axis=-1)/sigma**2)
    return np.concatenate([gp+gm, 1j*(gp-gm)], axis=-1)


def gaussian_pair_gram(centers,sigma):
    """Exact whole-space L2 Gram matrix for the real density coefficients.

    Parseval uses the exp(-2*pi*i*k.x) convention. This turns a coefficient
    sensitivity ball into a physically interpretable integrated density-energy
    ball. The imaginary coefficient of a center at zero is identically null.
    """
    centers=np.asarray(centers,dtype=float)
    sigma=np.broadcast_to(np.asarray(sigma,dtype=float),(len(centers),))
    variance=sigma[:,None]**2+sigma[None,:]**2
    scale=(2*np.pi*sigma[:,None]**2*sigma[None,:]**2/variance)**1.5
    same=scale*np.exp(-np.sum((centers[:,None]-centers[None,:])**2,axis=-1)/(2*variance))
    opposite=scale*np.exp(-np.sum((centers[:,None]+centers[None,:])**2,axis=-1)/(2*variance))
    zero=np.zeros_like(same)
    return np.block([[2*(same+opposite),zero],[zero,2*(same-opposite)]])


def density_energy_coordinates(centers,sigma,tolerance=1e-12):
    """Return T such that delta=T*b has density L2 energy ||b||^2.

    Numerically null eigenmodes are excluded from the declared model subspace.
    Record the cutoff and omitted modes; this is not a whole-space uncertainty
    guarantee about unrepresented densities.
    """
    gram=gaussian_pair_gram(centers,sigma)
    values,vectors=np.linalg.eigh(gram)
    keep=values>tolerance*values[-1]
    transform=vectors[:,keep]/np.sqrt(values[keep])[None,:]
    return transform,{'retained_modes':int(keep.sum()),'excluded_modes':int((~keep).sum()),
                      'relative_eigenvalue_cutoff':tolerance,'largest_eigenvalue':float(values[-1]),
                      'smallest_retained_eigenvalue':float(values[keep][0])}


def evaluate_pairs(k, centers, sigma, coefficient, gradient=False):
    centers=np.asarray(centers)
    sigma=np.broadcast_to(np.asarray(sigma), (len(centers),))
    c=np.asarray(coefficient[:len(centers)])+1j*np.asarray(coefficient[len(centers):])
    dp=np.asarray(k)[...,None,:]-centers
    dm=np.asarray(k)[...,None,:]+centers
    gp=np.exp(-0.5*np.sum(dp*dp,axis=-1)/sigma**2)
    gm=np.exp(-0.5*np.sum(dm*dm,axis=-1)/sigma**2)
    value=gp@c+gm@c.conjugate()
    if not gradient:
        return value
    grad=-np.sum(dp*(gp*c/sigma**2)[...,None]+dm*(gm*c.conjugate()/sigma**2)[...,None],axis=-2)
    return value,grad


def realify(z):
    return np.concatenate([np.real(z),np.imag(z)],axis=-1)


def image_design(k, ctf, centers, sigma, noise_std=1.0):
    f=gaussian_pair_features(k,centers,sigma)*np.asarray(ctf)[...,None]/noise_std
    return np.concatenate([f.real,f.imag],axis=-2)


def pose_jacobian(k, q, ctf, centers, sigma, pilot, box,
                  rotation_radius, shift_radius=0.0, noise_std=1.0):
    """Jacobian for right-multiplied row-vector rotations and image shifts.

    u in R^5 has ||u||<=1. Physical rotvec=rotation_radius*u[:3] (radians)
    and translation=shift_radius*u[3:] (working pixels).
    """
    f,g=evaluate_pairs(k,centers,sigma,pilot,gradient=True)
    rotations=np.cross(g,k)*rotation_radius
    translations=(-2j*np.pi/box)*np.asarray(q)*f[...,None]*shift_radius
    jac=np.concatenate([rotations,translations],axis=-1)*np.asarray(ctf)[...,None]/noise_std
    return np.concatenate([jac.real,jac.imag],axis=-2)


def perturbed_images(k, q, ctf, centers, sigma, coefficient, box, u,
                     rotation_radius, shift_radius=0.0, noise_std=1.0):
    rotations=Rotation.from_rotvec(np.asarray(u)[:,:3]*rotation_radius).as_matrix()
    kp=np.einsum('nqi,nij->nqj',k,rotations)
    phases=np.exp((-2j*np.pi/box)*np.einsum('nqi,ni->nq',q,np.asarray(u)[:,3:]*shift_radius))
    return evaluate_pairs(kp,centers,sigma,coefficient)*ctf*phases/noise_std


def pose_design_derivatives(k, q, ctf, centers, sigma, box,
                            rotation_radius, shift_radius=0.0, noise_std=1.0):
    """dA/du with shape (particles, realified pixels, coefficients, 5).

    For weights w, B*||sum_m w_m dA_m/du||_F bounds the first-order
    density-pose interaction for ||delta||<=B and ||u||<=1. A spectral norm
    would be tighter; the Frobenius relaxation permits group-norm optimization.
    """
    k,centers=np.asarray(k),np.asarray(centers)
    sigma=np.broadcast_to(np.asarray(sigma),(len(centers),))
    dp=k[...,None,:]-centers;dm=k[...,None,:]+centers
    gp=np.exp(-.5*np.sum(dp*dp,axis=-1)/sigma**2)
    gm=np.exp(-.5*np.sum(dm*dm,axis=-1)/sigma**2)
    gradp=-dp*(gp/sigma**2)[...,None]
    gradm=-dm*(gm/sigma**2)[...,None]
    gradient=np.concatenate([gradp+gradm,1j*(gradp-gradm)],axis=-2)
    features=np.concatenate([gp+gm,1j*(gp-gm)],axis=-1)
    rotations=np.cross(gradient,k[...,None,:])*rotation_radius
    translations=(-2j*np.pi/box)*np.asarray(q)[...,None,:]*features[...,None]*shift_radius
    derivatives=np.concatenate([rotations,translations],axis=-1)*np.asarray(ctf)[...,None,None]/noise_std
    return np.concatenate([derivatives.real,derivatives.imag],axis=-3)


def pose_remainder_bounds(k, q, ctf, centers, sigma, pilot, box,
                          rotation_radius, shift_radius, density_radius,
                          noise_std=1.0, structured=False):
    """Conservative analytic suprema over the entire normalized pose unit ball.

    Bound ||(A(u)-A(0))*delta|| by B*L and the pilot's second-order
    remainder by H/2. Returns particle-wise gamma=B*L+H/2, L and H.
    The full-Gaussian finite dictionary, fixed CTF, and common known real/imag
    noise standard deviation are part of the model; CTF errors are not included.
    With structured=True, retain the first-order density-pose interaction in
    the estimator and return gamma=(H+B*H_operator)/2, H_operator, H instead.
    """
    if not 0 <= rotation_radius <= np.pi or shift_radius < 0 or density_radius < 0:
        raise ValueError("Nonnegative radii and rotation radius <= pi required")
    centers=np.asarray(centers)
    sigma=np.broadcast_to(np.asarray(sigma), (len(centers),))
    c=np.asarray(pilot[:len(centers)])+1j*np.asarray(pilot[len(centers):])
    kr=np.linalg.norm(k,axis=-1)
    angular_displacement=2*kr*np.sin(rotation_radius/2)
    speed=rotation_radius*kr
    acceleration=rotation_radius**2*kr
    phase_speed=2*np.pi*np.linalg.norm(q,axis=-1)*shift_radius/box
    first=[]; second=[]
    for sign in [-1,1]:
        distance=np.linalg.norm(np.asarray(k)[...,None,:]+sign*centers,axis=-1)
        lo=np.maximum(0,distance-angular_displacement[...,None])
        hi=distance+angular_displacement[...,None]
        peak=np.clip(sigma,lo,hi)
        g0=np.exp(-0.5*(lo/sigma)**2)
        g1=peak/sigma**2*np.exp(-0.5*(peak/sigma)**2)
        # (1+r^2/sigma^2)*exp(-r^2/(2 sigma^2))/sigma^2 bounds
        # the Hessian operator norm and attains its maximum at r=sigma.
        g2=(1+(peak/sigma)**2)/sigma**2*np.exp(-0.5*(peak/sigma)**2)
        d1=speed[...,None]*g1+phase_speed[...,None]*g0
        d2=(speed[...,None]**2*g2+acceleration[...,None]*g1
            +2*(speed*phase_speed)[...,None]*g1+phase_speed[...,None]**2*g0)
        first.append(d1);second.append(d2)
    transfer=np.abs(ctf)/noise_std
    L=np.sqrt(np.sum(2*transfer**2*np.sum(first[0]**2+first[1]**2,axis=-1),axis=-1))
    h_per_frequency=transfer*np.sum((second[0]+second[1])*np.abs(c),axis=-1)
    H=np.linalg.norm(h_per_frequency,axis=-1)
    if structured:
        H_operator=np.sqrt(np.sum(2*transfer**2*np.sum(second[0]**2+second[1]**2,axis=-1),axis=-1))
        return .5*(H+density_radius*H_operator),H_operator,H
    return density_radius*L+0.5*H,L,H


def density_functionals(points, centers, sigma):
    """Continuous inverse-Fourier values; points use reciprocal frequency units."""
    points,centers=np.asarray(points),np.asarray(centers)
    sigma=np.broadcast_to(np.asarray(sigma),(len(centers),))
    phase=2*np.pi*(points@centers.T)
    envelope=2*(2*np.pi)**1.5*sigma**3*np.exp(-2*np.pi**2*np.sum(points*points,axis=-1)[:,None]*sigma**2)
    return np.concatenate([envelope*np.cos(phase),-envelope*np.sin(phase)],axis=-1)


def gaussian_average_functionals(points,centers,sigma,width):
    """Density averaged against normalized real-space Gaussian kernels.

    Points and width use the reciprocal units of the Fourier centers. For
    frequencies measured in field-of-view bins, multiply them by the physical
    field of view to get Angstroms. The integral is analytic, including the
    attenuation of each Fourier center; width=0 recovers point evaluation.
    """
    if width<0:raise ValueError('Nonnegative averaging width required')
    points=np.asarray(points);centers=np.asarray(centers)
    sigma=np.broadcast_to(np.asarray(sigma),(len(centers),))
    denominator=1+4*np.pi**2*width**2*sigma**2
    effective_sigma=sigma/np.sqrt(denominator)
    effective_centers=centers/denominator[:,None]
    attenuation=np.exp(-2*np.pi**2*width**2*np.sum(centers**2,axis=1)/denominator)
    values=density_functionals(points,effective_centers,effective_sigma)
    return values*np.tile(attenuation,2)
