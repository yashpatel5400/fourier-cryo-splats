"""Two-sided finite-orbit diagnostics for paired covariance tests.

No sampled orbit supremum is a continuous-pose certificate. All arithmetic
uses ordinary floating point, with explicit diagnostic rounding pads.
"""
import time
import numpy as np
from scipy.optimize import minimize, brentq
from scipy.stats import norm
from .uq_paired_covariance import direction_growth


def unrestricted_matrix_growth_upper(residual, variance_upper=1.):
    """Supremum of tr(T residual/v)+.5 logdet(I-T^2), -I<T<I.

    When residual=S-A and A is a nonnegative null covariance mixture, this
    is an upper bound for every feasible bilinear test, regardless of whether
    A minimizes distance to the cone. Spectral alignment gives the scalar
    formula g(d)=r-.5 log(1+r), r=(sqrt(1+4d^2)-1)/2.
    """
    delta=np.asarray(residual,float);v=float(variance_upper)
    if (delta.ndim!=2 or delta.shape[0]!=delta.shape[1] or not np.isfinite(delta).all()
            or not np.allclose(delta,delta.T,atol=1e-10,rtol=1e-10) or not np.isfinite(v) or v<=0):
        raise ValueError('Finite symmetric residual and positive variance bound required')
    eigenvalues=np.linalg.eigvalsh((delta+delta.T)/2)/v
    absolute=np.abs(eigenvalues)
    r=absolute*(2*absolute/(1+np.hypot(1,2*absolute)))
    raw=float(np.sum(r-.5*np.log1p(r)))
    pad=float(100*np.finfo(float).eps*len(delta)*max(1.,np.sum(absolute)))
    return dict(raw_spectral_value=raw,roundoff_diagnostic_pad=pad,upper=raw+pad)


def scalar_dual_tangent_upper(signal,penalty):
    """Correct root-evaluation shortfall with a full-domain concave tangent.

    For f concave on (-1,1), sup f<=f(x)+max(f'(x)(1-x),
    f'(x)(-1-x)). The root can be approximate; x=0 is a safe fallback.
    """
    s,a=float(signal),float(penalty)
    if not np.isfinite([s,a]).all() or min(s,a)<0:
        raise ValueError('Finite nonnegative inputs required')
    def derivative(x):return s-2*x/(1-x*x)-a/(1-x)**2
    try:x=float(brentq(derivative,-1+1e-14,1-1e-14,xtol=5e-15))
    except ValueError:x=0.
    terms=[s*x,np.log1p(-x*x),-a*x/(1-x)]
    value=float(sum(terms));d=derivative(x)
    tangent=max(d*(1-x),d*(-1-x),0.)
    pad=100*np.finfo(float).eps*(1+sum(abs(term) for term in terms)+abs(d))
    return dict(upper=max(0.,value+tangent)+pad,root=x,root_value=value,
        tangent_correction=float(tangent),roundoff_diagnostic_pad=float(pad))


def full_covariance_cone(means,signal_second_moment,maximum_seconds=120.,maximum_iterations=3000):
    """Kernel NNLS with retained feasible iterates and a repaired separator.

    ||sum lambda_i u_i u_i' - S||_F^2 uses kernel (u_i'u_j)^2.
    Nonnegative coefficients give a valid cone approximation even on timeout.
    The spectral upper and repaired feasible lower do not require optimizer
    convergence. Their gap measures the actual unresolved discrimination.
    """
    m=np.asarray(means,float);s=np.asarray(signal_second_moment,float)
    if m.ndim!=2 or s.shape!=(m.shape[1],)*2 or not np.isfinite(m).all() or not np.isfinite(s).all():
        raise ValueError('Finite matching means and signal covariance required')
    if not np.allclose(s,s.T,atol=1e-10) or np.linalg.eigvalsh(s).min()<-1e-9:
        raise ValueError('Positive-semidefinite signal second moment required')
    tick=time.perf_counter();norms=np.sum(m*m,axis=1);nonzero=norms>0
    u=m[nonzero]/np.sqrt(norms[nonzero,None]);scale=np.linalg.norm(s,'fro')
    if scale<=0 or not nonzero.any():raise ValueError('Nonzero signal and null required')
    target=s/scale;kernel=(u@u.T)**2
    linear=np.einsum('ni,ij,nj->n',u,target,u)
    initial_cov=m.T@m/len(m);amplitude=max(0.,float(np.sum(initial_cov*target)/np.sum(initial_cov**2)))
    initial=norms[nonzero]/len(m)*amplitude
    current=initial.copy();iterations=[0];history=[];best_value=[np.inf]
    def objective(x):
        product=kernel@x;value=.5*x@product-linear@x+.5
        return float(value),product-linear
    class TimeLimit(Exception):pass
    def callback(x):
        nonlocal current
        iterations[0]+=1;value,_=objective(x)
        if value<best_value[0]:current=np.maximum(x,0.).copy();best_value[0]=value
        if iterations[0]%100==0:history.append(dict(iteration=iterations[0],objective=value,seconds=time.perf_counter()-tick))
        if time.perf_counter()-tick>=maximum_seconds:raise TimeLimit
    try:
        fit=minimize(objective,initial,jac=True,method='L-BFGS-B',bounds=[(0,None)]*len(initial),
            callback=callback,options=dict(maxiter=int(maximum_iterations),gtol=1e-11,ftol=1e-15,maxls=50,maxcor=20))
        if objective(fit.x)[0]<=objective(current)[0]:current=np.maximum(fit.x,0.)
        status=int(fit.status);message=str(fit.message);success=bool(fit.success)
    except TimeLimit:
        status=99;message='Declared time limit; best feasible iterate retained';success=False
    coefficients=np.zeros(len(m));coefficients[nonzero]=current*scale/norms[nonzero]
    approximation=np.einsum('n,ni,nj->ij',coefficients,m,m,optimize=True)
    residual=(s-approximation);residual=(residual+residual.T)/2
    direction=residual/max(np.linalg.norm(residual,2),1e-30)
    violation=np.einsum('ni,ij,nj->n',u,direction,u,optimize=True)
    repair=max(0.,float(violation.max()))+1e-12
    direction-=repair*np.eye(len(s))
    max_repaired=float(np.max(np.einsum('ni,ij,nj->n',u,direction,u,optimize=True)))
    if max_repaired>1e-10:raise ArithmeticError('Numerical finite-view feasibility repair failed')
    values={};weights={}
    for v in [1.,2.,4.]:
        growth=direction_growth(direction,s,v);weights[str(v)]=growth.pop('weights')
        upper=unrestricted_matrix_growth_upper(residual,v)
        if growth['expected_log_lower']>upper['upper']+1e-8:raise ArithmeticError('Weak duality violated')
        values[str(v)]=dict(**growth,**upper,gap=upper['upper']-growth['expected_log_lower'])
    return dict(coefficients=coefficients,approximation=approximation,direction=direction,
        weights=weights,growth=values,residual_frobenius=float(np.linalg.norm(residual,'fro')),
        relative_residual_frobenius=float(np.linalg.norm(residual,'fro')/scale),
        status=status,message=message,solver_converged=success,iterations=iterations[0],history=history,
        mixture_mass=float(coefficients.sum()),active_views=int(np.sum(coefficients>1e-12)),
        maximum_normalized_raw_violation=float(violation.max()),direction_repair=repair,
        maximum_normalized_repaired_violation=max_repaired,
        seconds=time.perf_counter()-tick)


def log_growth_variance(means,weights,variance_upper=1.,sample_sizes=(1000,10000,100000)):
    """Exact first two log moments; normal-approximation power is only diagnostic.

    Both actual real-coordinate noise covariances are I, independent between
    exposures. weights acts on observations divided by sqrt(variance_upper).
    The viewing law is uniform on the supplied means. Repeated particles and
    two frames with unequal means are not covered.
    """
    m=np.asarray(means,float);t=np.asarray(weights,float);v=float(variance_upper)
    k=t/v;eigenvalues=np.linalg.eigvalsh(t)
    if np.max(abs(eigenvalues))>=1:raise ValueError('Moment domain violated')
    linear=m@k;conditional=np.sum(linear*m,axis=1)
    normalizer=.5*np.log1p(-eigenvalues**2).sum()
    mean=float(conditional.mean()+normalizer)
    variance=float(np.sum(k*k)+2*np.mean(np.sum(linear*linear,axis=1))+np.var(conditional))
    powers={str(n):float(norm.sf((np.log(20)-n*mean)/np.sqrt(n*variance))) if variance>0 else float(n*mean>=np.log(20)) for n in sample_sizes}
    return dict(mean=mean,variance=variance,normal_approximation_rejection_probability=powers,
        scope='Normal approximation at unit independent Gaussian exposure noise; not measured or guaranteed power.')
