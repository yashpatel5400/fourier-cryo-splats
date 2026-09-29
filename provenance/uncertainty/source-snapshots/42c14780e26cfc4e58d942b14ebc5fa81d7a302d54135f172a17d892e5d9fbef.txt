"""Target-directed reduced-space solves with full-space bias verification.

The certificate is about an explicitly declared ambient model, not merely the
initial Gaussian fitting dictionary. Enrichment uses only the operator and the
target; it never examines inference observations. This is a reduced-space convex
optimization construction; its general principle is not claimed to be novel.
"""
from dataclasses import dataclass
import numpy as np
from scipy.stats import norm
from scipy.sparse.linalg import LinearOperator,cg
from .uncertainty import optimize_certificate,bias_aware_half_width


def orthonormalize_columns(columns,tolerance=1e-11):
    """Rank-revealing orthonormal basis; omitted fitting modes remain in audit."""
    columns=np.asarray(columns,dtype=float)
    u,s,_=np.linalg.svd(columns,full_matrices=False)
    return u[:,s>tolerance*s[0]]


@dataclass
class SubspaceCertificate:
    weights:np.ndarray
    noise_sd:float
    ambient_bias:float
    restricted_bias:float
    objective:float
    dual_lower_bound:float
    half_width:float
    restricted_half_width:float
    basis_dimension:int
    converged:bool
    history:list


def enrich_certificate(operator,functional,initial_basis,density_radius,
                       alpha=.05,max_rounds=40,rtol=.005,inner_rtol=.0005,
                       inner_maxiter=250):
    """Solve a fixed-pose certificate while auditing a larger density space.

    S has orthonormal columns. Its restricted problem is a lower bound on the
    ambient optimum, while every w has a valid ambient objective using A* w.
    Add the component of ell-A* w orthogonal to S, with two reorthogonalizations.
    Store the best ambient feasible iterate and a monotone dual lower bound.
    A numerical NUFFT tolerance is not a formal interval-arithmetic certificate;
    reproducibility audits must report it and compare against direct sums.
    """
    ell=np.asarray(functional,dtype=float).ravel();B=float(density_radius)
    if ell.shape!=(operator.shape[1],) or B<=0 or max_rounds<1:
        raise ValueError('Compatible functional, positive radius and rounds required')
    S=orthonormalize_columns(initial_basis)
    if S.shape[0]!=len(ell) or not S.shape[1]:raise ValueError('A nonempty compatible initial basis is required')
    a=operator.forward_columns(S);z=norm.ppf(1-alpha/2)
    best_w=np.zeros(operator.shape[0]);best_cost=B*np.linalg.norm(ell);best_bias=best_cost
    best_restricted=best_cost;best_sd=0.;best_lower=0.;history=[];converged=False
    for step in range(max_rounds):
        reduced_ell=S.T@ell
        fit=optimize_certificate(a,np.zeros((1,len(a),0)),reduced_ell,B,maxiter=inner_maxiter,rtol=inner_rtol)
        residual=ell-operator.adjoint(fit.weights)
        bias=B*np.linalg.norm(residual);cost=z*fit.noise_sd+bias
        if cost<best_cost:
            best_w=fit.weights.copy();best_cost=cost;best_bias=bias;best_sd=fit.noise_sd
            best_restricted=B*np.linalg.norm(S.T@residual)
        # The restricted dual v lifts to S*v, preserving its norm and A*v.
        best_lower=max(best_lower,fit.dual_lower_bound)
        outside=residual-S@(S.T@residual)
        outside-=S@(S.T@outside)
        outside_norm=np.linalg.norm(outside)
        history.append({'round':step,'dimension':S.shape[1],'ambient_objective':float(cost),
                        'best_ambient_objective':float(best_cost),'dual_lower_bound':float(best_lower),
                        'relative_gap':float(max(0,best_cost-best_lower)/best_cost),
                        'restricted_objective':float(fit.objective),'outside_residual_norm':float(outside_norm),
                        'restricted_inner_converged':bool(fit.converged),'restricted_inner_iterations':fit.iterations,
                        'ambient_half_width':bias_aware_half_width(fit.noise_sd,bias,alpha),
                        'restricted_half_width':fit.half_width})
        converged=best_cost-best_lower<=rtol*best_cost
        if converged or outside_norm<=1e-12*max(np.linalg.norm(ell),1e-300):break
        if step+1<max_rounds:
            direction=outside/outside_norm
            S=np.column_stack([S,direction])
            a=np.column_stack([a,operator.forward(direction)])
    return SubspaceCertificate(best_w,best_sd,best_bias,best_restricted,best_cost,min(best_lower,best_cost),
        bias_aware_half_width(best_sd,best_bias,alpha),bias_aware_half_width(best_sd,best_restricted,alpha),
        S.shape[1],bool(converged),history)


def matrix_free_certificate(operator,functional,density_radius,alpha=.05,maxiter=100,
                            rtol=.005,cg_rtol=1e-7,cg_maxiter=500,smoothing=1e-7):
    """Ambient fixed-pose optimum using only forward/adjoint operations.

    Quadratic majorization gives (A*A+lambda I)u=ell and w=A u. For any u,
    v=u*min(B/||u||, z/||Au||) is dual feasible, even if CG has not converged.
    This is also a full-space comparator for Gaussian subspace enrichment.
    """
    ell=np.asarray(functional,dtype=float).ravel();B=float(density_radius)
    if ell.shape!=(operator.shape[1],) or B<=0 or maxiter<1:raise ValueError('Invalid functional, radius or iterations')
    z=norm.ppf(1-alpha/2);base=B*np.linalg.norm(ell);eps=max(base,1e-300)*smoothing
    w=np.zeros(operator.shape[0]);u=np.zeros(operator.shape[1]);residual=ell.copy()
    best_w=w.copy();best_cost=base;best_bias=base;best_sd=0.;best_lower=0.;history=[];converged=False
    for it in range(maxiter):
        sd=np.linalg.norm(w);bias=B*np.linalg.norm(residual)
        ridge=(z*z/np.hypot(z*sd,eps))/(B*B/np.hypot(bias,eps))
        op=LinearOperator((len(ell),len(ell)),matvec=lambda x:operator.adjoint(operator.forward(x))+ridge*x,dtype=float)
        u,info=cg(op,ell,x0=u,rtol=cg_rtol,atol=0,maxiter=cg_maxiter)
        w=operator.forward(u);residual=ell-operator.adjoint(w)
        sd=np.linalg.norm(w);bias=B*np.linalg.norm(residual);cost=z*sd+bias
        scale=min(B/max(np.linalg.norm(u),1e-300),z/max(sd,1e-300))
        lower=max(0.,float(scale*(ell@u)));best_lower=max(best_lower,lower)
        if cost<best_cost:best_cost=cost;best_w=w.copy();best_bias=bias;best_sd=sd
        history.append({'iteration':it+1,'ridge':float(ridge),'objective':float(cost),
                        'best_objective':float(best_cost),'lower_bound':float(best_lower),'cg_info':int(info),
                        'relative_gap':float(max(0,best_cost-best_lower)/max(best_cost,1e-300))})
        if best_cost-best_lower<=rtol*max(best_cost,1e-300):converged=True;break
    return SubspaceCertificate(best_w,float(best_sd),float(best_bias),float(best_bias),float(best_cost),
        min(best_lower,best_cost),bias_aware_half_width(best_sd,best_bias,alpha),
        bias_aware_half_width(best_sd,best_bias,alpha),operator.shape[1],converged,history)
