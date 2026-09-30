"""Constructive two-point ambiguity bounds for two allowed pose designs.

Classical Gaussian testing gives a fixed-length CI lower bound from each
feasible pair. Optimization improves the witness; it is not required for
validity. Continuous Gram integration pads are retained throughout the audit.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator,cg


class PhaseRotatedGram:
    """A translated Fourier observation operator is an orthogonal row rotation."""
    def __init__(self,gram,translations):
        self.gram=gram;self.phase=np.exp(-2j*np.pi*np.asarray(translations))
        if self.phase.shape!=(gram.n,gram.nq):raise ValueError('One phase per Fourier pair required')
        self.n,self.nq,self.shape=gram.n,gram.nq,gram.shape

    def rotate(self,vector,inverse=False):
        v=np.asarray(vector).reshape(self.n,2*self.nq)
        value=(v[:,:self.nq]+1j*v[:,self.nq:])*(self.phase.conj() if inverse else self.phase)
        return self.gram.pack(value.real,value.imag)

    def matvec(self,weights):return self.rotate(self.gram.matvec(self.rotate(weights,True)))
    def target(self,centers,signs,width):
        a,norm2=self.gram.target(centers,signs,width);return self.rotate(a),norm2
    def quadrature_error(self,weights):return self.gram.quadrature_error(self.rotate(weights,True))
    def preconditioner(self,ridge):
        original=self.gram.preconditioner(ridge)
        return LinearOperator(self.shape,matvec=lambda x:self.rotate(original@self.rotate(x,True)),dtype=float)


def two_pose_witness(grams,targets,target_norm2,pilot_means,weights,B,pilot_norm,tau):
    """Audit explicit continuous densities, with quadrature norm pads.

    Return both a feasible feature gap and a dual upper bound on the fixed-pair
    modulus. A data-space error pad includes errors in each Gram action.
    """
    if not np.isfinite([B,pilot_norm,tau,target_norm2]).all() or B<pilot_norm or pilot_norm<0 or tau<=0 or target_norm2<0:
        raise ValueError('Need B>=pilot norm>=0 and a positive finite distance budget')
    w=np.asarray(weights,float);residuals=[];actions=[];errors=[];roundoff_pads=[]
    if len(grams)!=2 or len(targets)!=2 or len(pilot_means)!=2:raise ValueError('Exactly two designs required')
    if not np.isfinite(w).all() or any(not np.isfinite(v).all() for v in [*targets,*pilot_means]):
        raise FloatingPointError('Nonfinite witness inputs')
    for gram,target in zip(grams,targets):
        action=gram.matvec(w);error=gram.quadrature_error(w)
        magnitude=max(1.,abs(target_norm2)+2*abs(w@target)+abs(w@action))
        rounding=64*np.finfo(float).eps*magnitude
        square=target_norm2-2*w@target+w@action+error['squared_field_norm']+rounding
        if not np.isfinite(square) or square<0:raise FloatingPointError('Invalid padded residual square')
        if any(not np.isfinite(v) or v<0 for v in error.values()):raise FloatingPointError('Invalid integration pad')
        residuals.append(np.sqrt(square))
        actions.append(action);errors.append(error);roundoff_pads.append(rounding)
    residuals=np.asarray(residuals);amplitude=np.divide(B,residuals,out=np.zeros_like(residuals),where=residuals>0)
    p=np.asarray(pilot_means[1])-np.asarray(pilot_means[0])
    mean=p+sum(a*(target-action) for a,target,action in zip(amplitude,targets,actions))
    pad=float(sum(a*error['gram_action_norm'] for a,error in zip(amplitude,errors)))
    distance_rounding=64*np.finfo(float).eps*max(1.,np.linalg.norm(p)+sum(a*(np.linalg.norm(t)+np.linalg.norm(g)) for a,t,g in zip(amplitude,targets,actions))+pad)
    distance=float(np.linalg.norm(mean)+pad+distance_rounding);scale=min(1.,tau/max(distance,np.finfo(float).tiny))
    feature=float(scale*sum(a*(target_norm2-w@target) for a,target in zip(amplitude,targets)))
    feature_pad=float(64*np.finfo(float).eps*max(1.,scale*sum(a*(abs(target_norm2)+np.sum(abs(w*t))) for a,t in zip(amplitude,targets))))
    upper_pad=float(64*np.finfo(float).eps*max(1.,B*residuals.sum()+np.sum(abs(w*p))+tau*np.linalg.norm(w)))
    upper=float(B*residuals.sum()-w@p+tau*np.linalg.norm(w)+upper_pad)
    if not np.isfinite([*amplitude,distance,scale,feature,upper,feature_pad]).all():raise FloatingPointError('Nonfinite witness result')
    gap=max(0.,abs(feature)-feature_pad)
    return {'feature_gap':gap,'signed_feature_gap':feature,'fixed_length_half_width_lower':gap/2,
            'fixed_pair_modulus_upper':upper,'mean_distance_upper':scale*distance,'unscaled_mean_distance_upper':distance,
            'common_density_scale':float(scale),'residual_amplitudes':amplitude.tolist(),
            'residual_norms_upper':residuals.tolist(),'mean_integration_pad':pad,
            'pilot_mean_difference_norm':float(np.linalg.norm(p)),'weight_norm':float(np.linalg.norm(w)),
            'feature_downward_roundoff_pad':feature_pad,'dual_upward_roundoff_pad':upper_pad,
            'mean_upward_roundoff_pad':float(distance_rounding),
            'residual_squared_roundoff_pads':roundoff_pads,'integration_errors':errors,
            'numerical_scope':'Analytic integration pads plus heuristic floating-point magnitude pads; not interval arithmetic or certified special-function rounding'}


def optimize_two_pose_modulus(grams,centers,signs,width,pilot_means,initial_weights,B,pilot_norm,tau,
                             maxiter=25,rtol=.005,callback=None,cg_maxiter=1000):
    """Iteratively reweighted quadratic minimization of the convex dual."""
    targets=[];norms=[]
    for gram in grams:
        a,norm2=gram.target(centers,signs,width);targets.append(a);norms.append(norm2)
    np.testing.assert_allclose(norms[0],norms[1]);target_norm2=float(norms[0])
    p=np.asarray(pilot_means[1])-np.asarray(pilot_means[0]);w=np.asarray(initial_weights,float).copy()
    best_lower=-np.inf;best_upper=np.inf;best_weights=None;history=[]
    for iteration in range(maxiter+1):
        witness=two_pose_witness(grams,targets,target_norm2,pilot_means,w,B,pilot_norm,tau)
        if witness['feature_gap']>best_lower:
            best_lower=witness['feature_gap'];best_weights=w.copy();best_witness=witness;lower_iteration=iteration
        if witness['fixed_pair_modulus_upper']<best_upper:
            best_upper=witness['fixed_pair_modulus_upper'];upper_weights=w.copy();upper_iteration=iteration
        if best_lower>best_upper+1e-7*max(1.,best_upper):raise AssertionError('Constructive lower exceeds dual upper')
        gap=max(0.,best_upper-best_lower)/max(best_upper,1e-30)
        row={'iteration':iteration,'dual_upper':witness['fixed_pair_modulus_upper'],
             'feasible_feature_gap':witness['feature_gap'],'best_lower':best_lower,'best_upper':best_upper,'relative_gap':gap}
        history.append(row)
        if callback is not None:callback(row)
        if gap<=rtol or iteration==maxiter:break
        alpha=B/np.maximum(witness['residual_norms_upper'],1e-12)
        ridge=tau/max(np.linalg.norm(w),1e-12)
        matrix=LinearOperator(grams[0].shape,matvec=lambda v:sum(a*g.matvec(v) for a,g in zip(alpha,grams))+ridge*v,dtype=float)
        base=grams[0].preconditioner(ridge/alpha.sum())
        preconditioner=LinearOperator(matrix.shape,matvec=lambda v:(base@v)/alpha.sum(),dtype=float)
        rhs=sum(a*t for a,t in zip(alpha,targets))+p
        count=[0]
        def tick(_):count[0]+=1
        w,info=cg(matrix,rhs,x0=w,M=preconditioner,rtol=1e-8,atol=0.,maxiter=cg_maxiter,callback=tick)
        row.update(cg_info=int(info),cg_iterations=count[0])
        if not np.isfinite(w).all():raise FloatingPointError('Nonfinite modulus iterate')
    return {'weights':best_weights,'upper_weights':upper_weights,'witness':best_witness,'dual_upper':best_upper,'feasible_lower':best_lower,
            'best_lower_iteration':lower_iteration,'best_upper_iteration':upper_iteration,
            'relative_modulus_gap':gap,'converged':bool(gap<=rtol),'history':history,'target_norm':np.sqrt(target_norm2),
            'distance_budget':tau,'scope':'Fixed-pair continuous modulus; lower bound remains valid without convergence. Selected poses need not maximize ambiguity.'}
