"""Candidate-only local perturbation scores, without a reference-map input.

The score design is a heuristic matched contrast, not a new estimator theorem.
Independent known-simulator calibration supplies its conditional testing claim.
"""
import numpy as np
from scipy.ndimage import gaussian_filter
from fourier_splats.uq_bispectrum import moment_features


def candidate_region(volume,field_A,width_A=20.):
    """Largest smoothed candidate-density peak in the central half-field cube."""
    rho=np.asarray(volume,float);box=rho.shape[0]
    if rho.shape!=(box,box,box) or not np.isfinite(rho).all() or field_A<=0 or width_A<=0:
        raise ValueError('Finite cubic candidate and positive physical scales required')
    grid=(np.arange(box)+.5)/box-.5;z,y,x=np.meshgrid(grid,grid,grid,indexing='ij')
    interior=(abs(x)<=.25)&(abs(y)<=.25)&(abs(z)<=.25)
    smooth=gaussian_filter(rho,width_A/field_A*box,mode='constant',cval=0.,truncate=4.)
    peak=np.unravel_index(np.argmax(np.where(interior,smooth,-np.inf)),rho.shape)
    center=np.array([grid[peak[2]],grid[peak[1]],grid[peak[0]]])
    mask=np.exp(-.5*((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)/(width_A/field_A)**2)
    return mask,dict(center_fraction_field=center.tolist(),center_index_zyx=list(map(int,peak)),
        width_A=width_A,field_A=field_A,smoothed_peak=float(smooth[peak]),
        removed_L2_fraction=float(np.linalg.norm(rho*mask)/np.linalg.norm(rho)))


def candidate_moment_direction(candidate_means,region_means,triads,remove_scale,design_removal=.25):
    """Match a prescribed 25% candidate-region deletion; optionally remove scale means."""
    m=np.asarray(candidate_means,complex);g=np.asarray(region_means,complex)
    if m.shape!=g.shape or m.ndim!=2:raise ValueError('Matching training-view Fourier means required')
    original=moment_features(m,triads).mean(axis=0)
    altered=moment_features(m-design_removal*g,triads).mean(axis=0)
    direction=altered-original;n=m.shape[1]
    if remove_scale:
        # Power and bispectrum means have disjoint coordinates and scale as a^2/a^3.
        for segment in [slice(0,n),slice(n,len(original))]:
            u=original[segment];energy=u@u
            if energy>0:direction[segment]-=(direction[segment]@u)/energy*u
    norm=np.linalg.norm(direction)
    if norm>1e-12:direction/=norm
    else:direction=np.zeros_like(direction)
    null_mean=float(original@direction);alternative_mean=float(altered@direction)
    return direction,dict(design_removal=design_removal,remove_scale=bool(remove_scale),
        original_direction_norm=float(norm),zero_direction=bool(not np.any(direction)),
        training_null_mean=null_mean,training_design_alternative_mean=alternative_mean,
        threshold=(null_mean+alternative_mean)/2,
        power_scale_mean=float(original[:n]@direction[:n]),bispectrum_scale_mean=float(original[n:]@direction[n:]))
