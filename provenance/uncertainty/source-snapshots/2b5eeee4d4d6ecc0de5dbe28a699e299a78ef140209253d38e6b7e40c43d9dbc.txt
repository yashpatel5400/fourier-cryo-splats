"""Envelope-mixture-guided refinement of a full SO(3) covering partition.

The refinement priority is heuristic. Validity comes from the complete cover,
cell upper kernels and classical scaled mixture dual, not that priority.
"""
import time
import numpy as np
from scipy.special import logsumexp
from .uq_continuous_mixture import FourierGaussianOrbit
from .uq_mixture_curvature import CurvatureGaussianOrbit
from .uq_mixture_fast import matrix_mixture_fit
from .uq_mixture_anchor import scaled_anchor_upper


class SelectiveCurvatureOrbit(CurvatureGaussianOrbit):
    """Avoid expensive curvature calculations outside a declared small box."""
    def cell(self, lower, upper):
        if np.sum(np.asarray(upper)-np.asarray(lower))/2 > np.deg2rad(5):
            return FourierGaussianOrbit.cell(self, lower, upper)
        return super().cell(lower, upper)


def refine_envelope_mixture(orbit, parent, *, max_splits=32768, batch_size=64,
                            mixture_iterations=100, lower_iterations=200,
                            lower_every=4, tolerance=1., wall_seconds=900.,
                            callback=None):
    """Continue a recorded covering partition, preserving replayable witnesses.

    ``parent`` holds arrays from continuous_mixture_upper. Its initial finite
    lower witness is retained; no oracle orientation or likelihood is used.
    Original log kernels, not floored matrix-EM quantities, certify bounds.
    """
    if (max_splits<0 or batch_size<1 or mixture_iterations<0 or lower_iterations<0
            or lower_every<1 or tolerance<=0 or wall_seconds<=0):
        raise ValueError('Invalid refinement budget')
    start=time.perf_counter()
    ids=np.array(parent['leaf_ids'],int).copy()
    low=np.array(parent['leaf_lower'],float).copy()
    high=np.array(parent['leaf_upper'],float).copy()
    upper=np.array(parent['leaf_log_envelopes'],float).copy()
    support=list(np.array(parent['support_log_kernels'],float))
    angles=list(np.array(parent['support_angles'],float))
    if (len(set(ids))!=len(ids) or low.shape!=(len(ids),3) or high.shape!=low.shape
            or upper.shape!=(len(ids),orbit.n) or np.any(high<low)
            or np.any(ids<0) or np.any(ids>=len(support))):
        raise ValueError('Incompatible initial cover and center arrays')
    count=int(parent['best_feasible_support_count'])
    feasible_ids=np.arange(count)
    feasible_weights=np.array(parent['best_feasible_weights'],float).copy()
    best_lower=float(logsumexp(np.array(support[:count]).T+
        np.log(feasible_weights)[None,:],axis=1).sum())
    best_upper=float('inf');best_cover=None;best_anchor=None
    weights=None;history=[];splits=[];total=0;batch=0;status='split_limit'
    while True:
        matrix=np.array(support)
        # At most four center candidates per image, plus the existing witness.
        # This restriction supplies only a feasible lower bound.
        if batch%lower_every==0 or total==max_splits:
            top=np.argpartition(matrix,-min(4,len(matrix)),axis=0)[-min(4,len(matrix)):]
            candidate=np.unique(top)
            finite=matrix_mixture_fit(matrix[candidate].T,
                max_iterations=lower_iterations,tolerance=.01)
            if finite['primal']>best_lower:
                best_lower=finite['primal'];feasible_ids=candidate
                feasible_weights=finite['weights'].copy()
        envelope=matrix_mixture_fit(upper.T,weights=weights,
            max_iterations=mixture_iterations,tolerance=.01)
        weights=envelope['weights'].copy()
        if envelope['upper']<best_upper:
            best_upper=envelope['upper'];best_anchor=envelope['upper_log_anchor'].copy()
            best_cover=dict(ids=ids.copy(),lower=low.copy(),upper=high.copy(),envelopes=upper.copy())
        if best_upper<best_lower-1e-7:
            raise FloatingPointError('Continuous upper below feasible lower')
        gap=max(best_upper-best_lower,0.)
        row=dict(splits=total,leaves=len(ids),support_count=len(support),
            lower=best_lower+orbit.normalization,upper=best_upper+orbit.normalization,
            gap=gap,envelope_fit_gap=envelope['gap'],
            envelope_fit_converged=envelope['converged'],
            envelope_fit_seconds=envelope['seconds'],seconds=time.perf_counter()-start)
        history.append(row)
        if callback is not None:callback(row)
        if gap<=tolerance:status='converged';break
        if total==max_splits:break
        if time.perf_counter()-start>=wall_seconds:status='wall_limit';break
        # Responsibility-weighted mismatch gives a mass-sensitive priority.
        # Residual dual violations provide a secondary exploration term.
        z=logsumexp(upper.T+np.log(weights)[None,:],axis=1)
        ratios=np.exp(np.minimum(upper-z[None,:],700.))
        slack=-np.expm1(np.minimum(matrix[ids]-upper,0.))
        mass_slack=weights*np.sum(ratios*slack,axis=1)
        score=ratios.sum(axis=1)
        priority=mass_slack+np.maximum(score-orbit.n,0.)/len(ids)+1e-12*score
        chosen=np.argsort(-priority,kind='stable')[:min(batch_size,max_splits-total,len(ids))]
        keep=np.ones(len(ids),bool);keep[chosen]=False
        child_ids=[];child_low=[];child_high=[];child_upper=[];child_weights=[]
        for position in chosen:
            axis=int(np.argmax(high[position]-low[position]))
            middle=(high[position,axis]+low[position,axis])/2
            left_high=high[position].copy();left_high[axis]=middle
            right_low=low[position].copy();right_low[axis]=middle
            children=[]
            for l,h in [(low[position],left_high),(right_low,high[position])]:
                cell=orbit.cell(l,h)
                # A parent upper is also an upper over either child.
                bound=np.minimum(cell['envelope'],upper[position])
                if np.any(bound<cell['exact']-1e-7):
                    raise FloatingPointError('Inherited envelope below child center')
                index=len(support);support.append(cell['exact']);angles.append(cell['middle'])
                children.append(index);child_ids.append(index);child_low.append(l.copy())
                child_high.append(h.copy());child_upper.append(bound)
                child_weights.append(weights[position]/2)
            splits.append([int(ids[position]),axis,*children])
        ids=np.concatenate([ids[keep],child_ids])
        low=np.concatenate([low[keep],child_low]);high=np.concatenate([high[keep],child_high])
        upper=np.concatenate([upper[keep],child_upper])
        weights=np.concatenate([weights[keep],child_weights]);weights/=weights.sum()
        total+=len(chosen);batch+=1
    # Retain the inherited monotone cover with any earlier valid anchor.
    replay=scaled_anchor_upper(upper.T,best_anchor)
    if replay<best_upper:
        best_upper=replay
        best_cover=dict(ids=ids.copy(),lower=low.copy(),upper=high.copy(),envelopes=upper.copy())
    if best_upper-best_lower<=tolerance:status='converged'
    return dict(log_likelihood_lower=best_lower+orbit.normalization,
        log_likelihood_upper=best_upper+orbit.normalization,
        gap=max(best_upper-best_lower,0.),converged=status=='converged',status=status,
        seconds=time.perf_counter()-start,history=history,
        split_history=np.asarray(splits,int).reshape(-1,4),
        leaf_ids=ids,leaf_lower=low,leaf_upper=high,leaf_log_envelopes=upper,
        support_angles=np.array(angles),support_log_kernels=np.array(support),
        best_anchor_log_z=best_anchor,best_cover_ids=best_cover['ids'],
        best_cover_lower=best_cover['lower'],best_cover_upper=best_cover['upper'],
        best_cover_log_envelopes=best_cover['envelopes'],
        best_feasible_ids=feasible_ids,best_feasible_weights=feasible_weights,
        normalization=orbit.normalization,
        scope='Continuous mixture computation, supplied noise/CTFs, zero shifts; no predictive-test outcome.')
