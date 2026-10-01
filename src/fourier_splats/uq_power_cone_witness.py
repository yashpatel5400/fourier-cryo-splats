"""Nonnegative finite-orientation witnesses for linear cross-power limitations."""
import numpy as np
from scipy.optimize import linprog
from .uq_power_cone import scalar_dual_maximum


def cone_witness(powers, signal):
    p, s = np.asarray(powers,float), np.asarray(signal,float)
    if (p.ndim!=2 or s.shape!=(p.shape[1],) or not np.isfinite(p).all()
            or not np.isfinite(s).all() or min(p.min(),s.min())<0 or s.max()<=0):
        raise ValueError('Nonnegative finite powers and nonzero target required')
    scale=np.maximum(s,1e-6*s.max());matrix=p.T/scale[:,None];target=s/scale
    options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9,time_limit=60.)
    first=linprog(np.zeros(len(p)),A_eq=matrix,b_eq=target,bounds=(0,None),method='highs',options=options)
    second=None
    if first.success:
        coefficients=np.maximum(first.x,0.)
    else:
        inequalities=np.vstack([np.column_stack([matrix,-np.ones(len(s))]),
                                 np.column_stack([-matrix,-np.ones(len(s))])])
        objective=np.zeros(len(p)+1);objective[-1]=1.
        second=linprog(objective,A_ub=inequalities,b_ub=np.concatenate([target,-target]),
            bounds=(0,None),method='highs',options=options)
        if not second.success:
            raise ArithmeticError(f'Equality status {first.status}; approximation status {second.status}')
        coefficients=np.maximum(second.x[:-1],0.)
    approximation=p.T@coefficients;residual=approximation-s
    return dict(coefficients=coefficients,approximation=approximation,
        equality_status=int(first.status),equality_success=bool(first.success),equality_message=first.message,
        approximation_status=None if second is None else int(second.status),
        maximum_scaled_residual=float(np.max(abs(residual)/scale)),
        absolute_residual_norm=float(np.linalg.norm(residual)),
        mixture_mass=float(coefficients.sum()),active_views=int(np.sum(coefficients>1e-12)))


def power_witness_log_upper(signal, approximation, transfer_squared, variance_upper=1.):
    """Numerical separable dual using a nonnegative cone approximation.

    Caller must establish approximation = P'lambda with lambda>=0 for null
    candidate powers P. The source map/operator is not checked here.
    """
    s,a,c=np.broadcast_arrays(np.asarray(signal,float),np.asarray(approximation,float),
        np.asarray(transfer_squared,float))
    v=float(variance_upper)
    if not all(np.isfinite(x).all() and np.all(x>=0) for x in [s,a,c]) or not np.isfinite(v) or v<=0:
        raise ValueError('Finite nonnegative inputs and positive variance required')
    values=[scalar_dual_maximum(x,y)[0] for x,y in zip((s*c/v).ravel(),(a*c/v).ravel())]
    return float(sum(values))
