"""Classical folded-normal width search on a fixed-design ridge path.

Analytic interval bounds use monotonicity of the exact ridge tradeoff and a
CG residual bound. Floating point and NUFFT remain diagnostically padded,
not validated interval arithmetic. No pose-estimation guarantee is added.
"""
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg
from .uq_intervals import bias_aware_half_width_stable


def folded_width_derivatives(sd, bias, alpha=.05):
    if sd <= 0 or bias < 0 or not 0 < alpha < .5:
        raise ValueError('Positive SD, nonnegative bias, and alpha in (0,.5) required')
    q=bias_aware_half_width_stable(sd,bias,alpha)
    derivative_bias=np.tanh(q*bias/(sd*sd))
    derivative_sd=(q-bias*derivative_bias)/sd
    return float(derivative_sd),float(derivative_bias)


def folded_ridge_search(gram, centers, signs, width, density_radius, *, alpha=.05,
        initial_ridge=1., relative_tolerance=.005, maximum_evaluations=80,
        maximum_seconds=180., cg_rtol=1e-10, cg_maxiter=2000):
    a,ell2=gram.target(centers,signs,width);a=np.asarray(a,float)
    B=float(density_radius);initial_ridge=float(initial_ridge)
    if (B<=0 or ell2<=0 or not np.isfinite([B,ell2,initial_ridge]).all() or initial_ridge<=0
            or not 0<alpha<.5 or relative_tolerance<=0 or maximum_evaluations<3):
        raise ValueError('Positive finite scales, valid alpha, and at least three evaluations required')
    start=time.perf_counter();cache={};history=[]
    base=B*np.sqrt(ell2)
    best=dict(weights=np.zeros_like(a),half_width=float(base),noise_sd=0.,bias=float(base),ridge=None)
    quadrature=getattr(gram,'quadrature_error',lambda w:{'squared_field_norm':0.,'gram_action_norm':0.})

    def evaluate(log_ridge):
        nonlocal best
        key=float(log_ridge)
        if key in cache:return cache[key]
        ridge=float(np.exp(key));system=LinearOperator(gram.shape,matvec=lambda w:gram.matvec(w)+ridge*w,dtype=float)
        pre=gram.preconditioner(ridge) if hasattr(gram,'preconditioner') else LinearOperator(
            gram.shape,matvec=lambda w:w/(gram.diagonal+ridge),dtype=float)
        warm=None if not cache else cache[min(cache,key=lambda old:abs(old-key))]['weights']
        iterations=[0]
        def count(_):iterations[0]+=1
        w,info=cg(system,a,x0=warm,M=pre,rtol=cg_rtol,atol=0,maxiter=cg_maxiter,callback=count)
        gw=gram.matvec(w);errors=quadrature(w)
        sd=float(np.linalg.norm(w));h2=float(ell2-2*w@a+w@gw)
        pad=float(50*np.finfo(float).eps*len(w)*(ell2+2*abs(w@a)+abs(w@gw)))
        norm_error=float(errors['squared_field_norm'])+pad
        if h2 < -norm_error:raise ArithmeticError('Negative residual norm beyond diagnostic pad')
        h_upper=np.sqrt(max(0.,h2+norm_error));h_lower=np.sqrt(max(0.,h2-norm_error))
        solve_residual=float(np.linalg.norm(a-gw-ridge*w))+float(errors['gram_action_norm'])
        # For G>=0: ||(G+lambda I)^-1||<=1/lambda, and
        # ||G^(1/2)(G+lambda I)^-1||<=1/(2 sqrt(lambda)).
        error_w=solve_residual/ridge;error_h=solve_residual/(2*np.sqrt(ridge))
        sd_lower=max(0.,sd-error_w);bias_lower=B*max(0.,h_lower-error_h)
        bias=B*h_upper;half=bias_aware_half_width_stable(sd,bias,alpha)
        row=dict(log_ridge=key,ridge=ridge,half_width=half,noise_sd=sd,bias=float(bias),
            exact_path_noise_sd_lower=sd_lower,exact_path_bias_lower=bias_lower,
            cg_info=int(info),cg_iterations=iterations[0],solve_residual_upper=solve_residual,
            weight_error_upper=error_w,residual_function_error_upper=error_h,
            norm2_roundoff_diagnostic_pad=pad,quadrature_errors=errors,
            elapsed_seconds=time.perf_counter()-start)
        cache[key]=dict(weights=w,**row);history.append(row)
        if half<best['half_width']:
            best=dict(weights=w.copy(),half_width=half,noise_sd=sd,bias=float(bias),ridge=ridge)
        return cache[key]

    def lower_intervals():
        keys=sorted(cache);rows=[]
        # q(s,b) is coordinatewise increasing for alpha<1/2. Along the
        # exact path, s decreases and b increases with lambda.
        first,last=cache[keys[0]],cache[keys[-1]]
        rows.append(dict(left=None,right=keys[0],lower=bias_aware_half_width_stable(
            first['exact_path_noise_sd_lower'],0.,alpha)))
        for left,right in zip(keys[:-1],keys[1:]):
            rows.append(dict(left=left,right=right,lower=bias_aware_half_width_stable(
                cache[right]['exact_path_noise_sd_lower'],cache[left]['exact_path_bias_lower'],alpha)))
        rows.append(dict(left=keys[-1],right=None,lower=last['exact_path_bias_lower']))
        # This numerical downward pad is not a validated arithmetic enclosure.
        for row in rows:row['lower']=max(0.,row['lower']-100*np.finfo(float).eps*max(base,1.))
        return rows

    # The declared finite range is searched to tolerance. The two outer-ray
    # bounds are retained separately: an unobserved target component can make
    # them too loose to prove global optimality, even in a small dense model.
    origin=np.log(initial_ridge)
    for offset in [-np.log(64.),0.,np.log(64.)]:evaluate(origin+offset)
    while True:
        intervals=lower_intervals();worst=min(intervals[1:-1],key=lambda r:r['lower'])
        gap=max(0.,best['half_width']-worst['lower'])/best['half_width']
        if gap<=relative_tolerance:reason='relative_gap';break
        if len(cache)>=maximum_evaluations:reason='evaluation_limit';break
        if time.perf_counter()-start>=maximum_seconds:reason='time_limit';break
        left,right=worst['left'],worst['right'];key=(left+right)/2
        evaluate(key)
    global_lower=float(min(base,min(r['lower'] for r in intervals)))
    best.update(lower_bound=float(min(base,min(r['lower'] for r in intervals[1:-1]))),relative_gap=float(gap),
        global_lower_bound=global_lower,
        global_relative_gap=float(max(0.,best['half_width']-global_lower)/best['half_width']),
        searched_ridge_interval=[initial_ridge/64,initial_ridge*64],
        converged=bool(gap<=relative_tolerance),stop_reason=reason,history=history,
        intervals=intervals,evaluations=len(cache),seconds=time.perf_counter()-start,
        no_data_half_width=float(base),all_cg_converged=all(r['cg_info']==0 for r in history),
        scope='Fixed-design folded-width search on the declared finite ridge interval; global outer-ray bounds reported separately. Diagnostic floating-point padding.')
    return best
