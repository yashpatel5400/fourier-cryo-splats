"""Oracle finite-view power-cone optimization, not a continuous pose certificate."""
import numpy as np
import cvxpy as cp
from scipy.optimize import brentq


def scalar_dual_maximum(signal, penalty):
    """Supremum on -1<x<1 of sx+log(1-x²)-a*x/(1-x), a>=0.

    Strict concavity gives a unique derivative root. This is floating-point
    evaluation, not validated arithmetic. The box used for the primal solve
    does not restrict this dual maximization.
    """
    s, a = float(signal), float(penalty)
    if not np.isfinite([s, a]).all() or min(s, a) < 0:
        raise ValueError('Finite nonnegative signal and multiplier penalty required')
    def derivative(x):
        return s-2*x/(1-x*x)-a/(1-x)**2
    x = brentq(derivative, -1+1e-14, 1-1e-14, xtol=5e-15)
    value = s*x+np.log1p(-x*x)-a*x/(1-x)
    # x=0 is always feasible and has value zero.
    return float(max(value, 0.)), float(x)


class FinitePowerCone:
    def __init__(self, views, frequencies):
        self.x = cp.Variable(frequencies)
        self.p = cp.Parameter((views, frequencies), nonneg=True)
        self.signal = cp.Parameter(frequencies, nonneg=True)
        self.constraint = self.p@(cp.inv_pos(1-self.x)-1) <= 0
        objective = self.signal@self.x+cp.sum(cp.log(1-cp.square(self.x)))
        self.problem = cp.Problem(cp.Maximize(objective),
            [self.constraint, self.x <= .95, self.x >= -.95])

    def solve(self, powers, signal_power, variance_upper=1.):
        p, s = np.asarray(powers,float), np.asarray(signal_power,float)
        v=float(variance_upper)
        if (p.shape!=self.p.shape or s.shape!=self.signal.shape
                or not np.isfinite(v) or v<=0 or not np.isfinite(p).all()
                or not np.isfinite(s).all() or min(p.min(),s.min())<0):
            raise ValueError('Matching nonnegative powers and a positive variance bound required')
        totals=p.sum(axis=1)
        normalized=p/np.where(totals>0,totals,1.)[:,None]
        self.p.value=normalized;self.signal.value=s/v
        self.problem.solve(solver='CLARABEL', max_iter=300, tol_gap_abs=1e-8,
            tol_gap_rel=1e-8,tol_feas=1e-8,warm_start=False)
        if self.x.value is None or self.constraint.dual_value is None:
            raise ArithmeticError(f'No primal/dual iterate: {self.problem.status}')
        raw=np.asarray(self.x.value,float);multipliers=np.maximum(self.constraint.dual_value,0.)
        if not np.isfinite(raw).all() or np.max(abs(raw))>=1:
            raise ArithmeticError('Primal outside Gaussian moment domain')
        # Shrinking makes room inside the solve box before a downward repair
        # of the linearized coefficients. Subtracting the row-normalized
        # violation ensures every nonzero normalized row is feasible.
        x=raw*(1-1e-6);u=x/(1-x)
        correction=max(0.,float((normalized@u).max()))+1e-12
        u-=correction;x=u/(1+u)
        fallback=not np.isfinite(x).all() or np.max(abs(x))>=1
        if fallback:x=np.zeros_like(raw)
        residual=float((normalized@(x/(1-x))).max())
        if residual>1e-10:raise ArithmeticError('Power-cone repair failed')
        dual_terms=[scalar_dual_maximum(a,b)[0] for a,b in zip(s/v,normalized.T@multipliers)]
        upper=float(np.sum(dual_terms))
        objective=lambda w:float(s@w/v+np.log1p(-w*w).sum())
        lower=objective(x)
        # Zero is always a feasible candidate; preserve the raw iterate even
        # when the repaired solution is inferior by numerical tolerance.
        selected=np.zeros_like(x) if lower<0 else x
        lower=max(0.,lower)
        if upper<lower-1e-7*max(1.,abs(lower)):
            raise ArithmeticError('Numerical weak-duality violation')
        return dict(weights=selected/v, raw_dimensionless_weights=raw,
            multipliers=multipliers, expected_log_lower=lower,
            expected_log_dual_upper=upper, gap=max(0.,upper-lower),
            raw_expected_log=objective(raw),status=self.problem.status,
            iterations=self.problem.solver_stats.num_iters,
            maximum_normalized_violation=residual,
            coefficient_repair=correction,repair_fallback=fallback,
            maximum_raw_abs_weight=float(np.max(abs(raw))),
            full_domain_dual=True,variance_upper=v)
