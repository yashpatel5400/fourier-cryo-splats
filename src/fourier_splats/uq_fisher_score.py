"""Classical covariance-aware, scale-constrained moment discriminants.

This is a finite-training Fisher/Rayleigh construction. Calibration, not a
Gaussian feature approximation, determines the independent test guarantee.
"""
import numpy as np
from scipy.linalg import cho_factor,cho_solve


def constrained_fisher_direction(null_mean,alternative_mean,covariance,power_dimensions,shrinkage=.1):
    """Optimize squared mean gap / regularized variance subject to scale constraints.

S = (1-shrinkage)*covariance + shrinkage*trace(covariance)/p*I.
The power and bispectrum mean blocks separately have zero score weight.
The unit-Euclidean direction is oriented toward the declared alternative.
"""
    u=np.asarray(null_mean,float);v=np.asarray(alternative_mean,float);s=np.asarray(covariance,float)
    p=len(u);n=int(power_dimensions)
    if u.ndim!=1 or v.shape!=u.shape or s.shape!=(p,p) or not 0<n<=p or not 0<shrinkage<=1:
        raise ValueError('Matching moments, square covariance, block size and positive shrinkage required')
    if not all(np.isfinite(x).all() for x in [u,v,s]) or not np.allclose(s,s.T,atol=1e-10,rtol=1e-10):
        raise ValueError('Finite means and symmetric covariance required')
    average=float(np.trace(s)/p)
    if average<=0:raise ValueError('Positive average covariance eigenvalue required')
    reg=(1-shrinkage)*s+shrinkage*average*np.eye(p)
    factor=cho_factor(reg,lower=True)
    blocks=[]
    for sl in [slice(0,n),slice(n,p)]:
        c=np.zeros(p);c[sl]=u[sl];norm=np.linalg.norm(c)
        if norm>0:blocks.append(c/norm)
    constraints=np.column_stack(blocks) if blocks else np.empty((p,0))
    gap=v-u;invgap=cho_solve(factor,gap)
    if len(blocks):
        invu=cho_solve(factor,constraints)
        dual=np.linalg.solve(constraints.T@invu,constraints.T@invgap)
        raw=invgap-invu@dual
    else:raw=invgap
    norm=float(np.linalg.norm(raw))
    direction=raw/norm if norm>1e-12 else np.zeros_like(raw)
    null_score=float(u@direction);alt_score=float(v@direction)
    variance=float(direction@reg@direction)
    return direction,dict(shrinkage=float(shrinkage),covariance_trace=float(np.trace(s)),
        original_direction_norm=norm,zero_direction=bool(not np.any(direction)),
        training_null_mean=null_score,training_design_alternative_mean=alt_score,
        threshold=(null_score+alt_score)/2,power_scale_mean=float(u[:n]@direction[:n]),
        bispectrum_scale_mean=float(u[n:]@direction[n:]),regularized_training_variance=variance,
        regularized_squared_separation=(alt_score-null_score)**2/variance if variance>0 else 0.,
        maximum_constraint_residual=float(np.max(abs(constraints.T@direction))) if len(blocks) else 0.)
