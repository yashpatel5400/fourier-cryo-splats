"""Matched linear-Gaussian reference intervals with explicit statistical targets.

These are transparent mathematical baselines, not executions of external
cryo-EM programs. A posterior credible interval promises prior-predictive
calibration under its specified prior, not uniform coverage over a fixed ball.
"""
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.stats import norm
from scipy.sparse.linalg import LinearOperator,cg,aslinearoperator


def gaussian_reference(a,ell,prior_sd,j=None,pose_sd=0.,alpha=.05):
    a=np.asarray(a);ell=np.asarray(ell)
    if prior_sd<=0 or pose_sd<0:raise ValueError('Positive density prior scale required')
    if j is None or pose_sd==0:
        vinv_a=a
    else:
        n,m,q=j.shape;ab=a.reshape(n,m,-1)
        inv=np.linalg.inv(np.eye(q)[None]+pose_sd**2*j.transpose(0,2,1)@j)
        vinv_a=(ab-pose_sd**2*j@inv@j.transpose(0,2,1)@ab).reshape(a.shape)
    precision=a.T@vinv_a+np.eye(a.shape[1])/prior_sd**2
    factor=cho_factor(precision)
    v=cho_solve(factor,ell);w=vinv_a@v
    z=norm.ppf(1-alpha/2)
    return {'weights':w,'posterior_half_width':float(z*np.sqrt(ell@v)),
            'diagonal_vi_half_width':float(z*np.sqrt(np.sum(ell**2/np.diag(precision)))),
            'measurement_noise_half_width':float(z*np.linalg.norm(w)),
            'density_bias_direction':a.T@w-ell}


def gaussian_reference_operator(a,ell,prior_sd,j=None,pose_sd=0.,alpha=.05,
                                gram_diagonal=None,rtol=1e-9,maxiter=2000):
    """Full-covariance linear Gaussian reference without assembling A or A'A.

    A scalar isotropic density prior is explicit. Optional Gaussian pose error
    enters the additive linearized model; nonlinear interactions are not covered.
    """
    a=aslinearoperator(a);ell=np.asarray(ell,dtype=float)
    if prior_sd<=0 or pose_sd<0:raise ValueError('Positive density prior scale required')
    if j is None or pose_sd==0:
        def inverse_noise(v):return v
    else:
        n,m,q=j.shape
        inverse=np.linalg.inv(np.eye(q)[None]+pose_sd**2*j.transpose(0,2,1)@j)
        def inverse_noise(v):
            v=v.reshape(n,m)
            return (v-pose_sd**2*np.einsum('nmq,nqr,nsr,ns->nm',j,inverse,j,v)).ravel()
    ridge=1/prior_sd**2
    precision=LinearOperator((a.shape[1],a.shape[1]),matvec=lambda v:a.rmatvec(inverse_noise(a.matvec(v)))+ridge*v,dtype=float)
    precondition=None
    if gram_diagonal is not None:
        diagonal=np.broadcast_to(gram_diagonal,ell.shape)+ridge
        precondition=LinearOperator(precision.shape,matvec=lambda v:v/diagonal,dtype=float)
    v,info=cg(precision,ell,M=precondition,rtol=rtol,atol=0,maxiter=maxiter)
    residual=float(np.linalg.norm(precision@v-ell)/max(np.linalg.norm(ell),1e-300))
    if info or residual>max(10*rtol,1e-8):raise RuntimeError(f'Gaussian reference CG did not converge: {info}, {residual}')
    w=inverse_noise(a@v);variance=float(ell@v)
    if variance<0:raise ArithmeticError('Negative posterior variance')
    z=norm.ppf(1-alpha/2)
    return {'weights':w,'posterior_half_width':float(z*np.sqrt(variance)),
            'measurement_noise_half_width':float(z*np.linalg.norm(w)),
            'density_bias_direction':a.rmatvec(w)-ell,'relative_normal_residual':residual}


def prior_predictive_coverage(a,ell,prior_sd,reference,j=None,pose_sd=0.,alpha=.05):
    """Exact integrated coverage under the same linear Gaussian prior/model."""
    w=reference['weights'];bias_direction=a.T@w-ell
    variance=float(prior_sd**2*(bias_direction@bias_direction)+w@w)
    if j is not None and pose_sd:
        jt=np.einsum('nmq,nm->nq',j,w.reshape(j.shape[:2]))
        variance+=pose_sd**2*float(np.sum(jt*jt))
    half=reference['posterior_half_width']
    return float(2*norm.cdf(half/np.sqrt(variance))-1)


def unregularized_reference(a,ell,rtol=1e-10,alpha=.05):
    """Sampling interval only when the requested functional is identifiable."""
    u,s,vh=np.linalg.svd(a,full_matrices=False)
    keep=s>rtol*s[0]
    projected=vh[keep].T@(vh[keep]@ell)
    defect=float(np.linalg.norm(ell-projected))
    w=u[:,keep]@((vh[keep]@ell)/s[keep])
    identifiable=defect<=rtol*max(1.,np.linalg.norm(ell))
    return {'weights':w,'half_width':float(norm.ppf(1-alpha/2)*np.linalg.norm(w)) if identifiable else None,
            'identifiable':bool(identifiable),'functional_defect':defect,'rank':int(keep.sum()),
            'smallest_retained_singular_value':float(s[keep][-1])}
