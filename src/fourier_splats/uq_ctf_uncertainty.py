"""Uniform CTF/gain/envelope sensitivities for continuous density audits.

The declared radii are assumptions, not fitted confidence regions. Voltage,
spherical aberration and amplitude contrast stay fixed. Defocus/astigmatism,
phase, relative gain and relative B-factor uncertainties are bounded below.
"""
import numpy as np


def ctf_uncertainty_envelope(q_inverse_A,parameters,defocus_radius_A=0.,
                             astigmatism_angle_radius_degrees=0.,phase_radius_degrees=0.,
                             relative_gain_radius=0.,relative_B_radius_A2=0.):
    q=np.asarray(q_inverse_A,dtype=float);p=np.asarray(parameters,dtype=float)
    if q.ndim!=2 or q.shape[1]!=2 or p.ndim!=2 or p.shape[1]!=9:
        raise ValueError('Expected q by 2 frequencies and particle by 9 CTF parameters')
    if not np.isfinite(q).all() or not np.isfinite(p).all() or (p[:,5]<=0).any() or (abs(p[:,7])>1).any():
        raise ValueError('Finite physical CTF parameters required')
    n=len(p)
    d,a,ph,g,b=[np.broadcast_to(x,(n,)).astype(float)[:,None] for x in
                 [defocus_radius_A,astigmatism_angle_radius_degrees,phase_radius_degrees,
                  relative_gain_radius,relative_B_radius_A2]]
    if any(not np.isfinite(x).all() or (x<0).any() for x in [d,a,ph,g,b]) or (g>=1).any():
        raise ValueError('Nonnegative finite radii and gain radius below one required')
    s2=np.sum(q*q,axis=1)[None];az=np.arctan2(q[:,1],q[:,0])[None]
    voltage=p[:,5,None]*1000;lam=12.2639/np.sqrt(voltage+.97845e-6*voltage**2)
    du,dv,angle,cs,ac,phase=[p[:,j,None] for j in [2,3,4,6,7,8]]
    df=.5*(du+dv+(du-dv)*np.cos(2*(az-np.deg2rad(angle))))
    gamma=-np.pi*df*lam*s2+.5*np.pi*cs*1e7*lam**3*s2*s2-np.deg2rad(phase)
    theta=gamma-np.arcsin(ac);nominal=np.sin(theta)
    # Vary defocus axes at the perturbed angle, then vary the nominal angle.
    # The resulting bound also covers simultaneous changes of both axes/angle.
    df_radius=d+abs(du-dv)*np.sin(np.minimum(np.deg2rad(a),np.pi/2))
    phase_radius=np.pi*lam*s2*df_radius+np.deg2rad(ph)
    lo,hi=theta-phase_radius,theta+phase_radius
    minimum=np.minimum(np.sin(lo),np.sin(hi));maximum=np.maximum(np.sin(lo),np.sin(hi))
    has_max=np.ceil((lo-np.pi/2)/(2*np.pi))<=np.floor((hi-np.pi/2)/(2*np.pi))
    has_min=np.ceil((lo+np.pi/2)/(2*np.pi))<=np.floor((hi+np.pi/2)/(2*np.pi))
    maximum=np.where(has_max,1.,maximum);minimum=np.where(has_min,-1.,minimum)
    lower_amplitude=(1-g)*np.exp(-b*s2/4);upper_amplitude=(1+g)*np.exp(b*s2/4)
    products=np.stack([lower_amplitude*minimum,lower_amplitude*maximum,
                       upper_amplitude*minimum,upper_amplitude*maximum])
    transfer_min=products.min(axis=0);transfer_max=products.max(axis=0)
    envelope=np.maximum(abs(transfer_min-nominal),abs(transfer_max-nominal))
    return {'absolute_transfer_error':envelope,'nominal_transfer_float64':nominal,
            'minimum_transfer':transfer_min,'maximum_transfer':transfer_max,
            'phase_radius_radians':phase_radius}


def continuous_ctf_bias_bound(weights,noise_std,error_envelope,total_density_norm):
    """Triangle/Cauchy bound valid jointly with arbitrary allowed poses.

    On the unit-volume cube, ||Re(c exp(i phase))||_2 <= |c| for every pose.
    Add this term to a nominal-CTF density/pose bound. No Taylor truncation or
    division by a CTF zero occurs. This intentionally loses cross-frequency
    cancellation; it may be very conservative.
    """
    envelope=np.asarray(error_envelope,dtype=float)
    if envelope.ndim!=2 or not np.isfinite(envelope).all() or (envelope<0).any():
        raise ValueError('Finite nonnegative transfer error matrix required')
    if not np.isfinite(noise_std) or noise_std<=0 or not np.isfinite(total_density_norm) or total_density_norm<0:
        raise ValueError('Positive noise and nonnegative finite density norm required')
    n,nq=envelope.shape;w=np.asarray(weights,dtype=float).reshape(n,2*nq)
    return float(total_density_norm/noise_std*np.sum(envelope*np.hypot(w[:,:nq],w[:,nq:])))
