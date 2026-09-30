"""Small conic subproblems for spectral support exchange.

Every cut is supplied as a full-weight linear support of the quadrature pose
norm. The restricted model is an optimization device. Final continuous primal
and full-space dual audits are separate; a restricted conic optimum alone is
not an interval certificate or a full-space optimum.
"""
import numpy as np
import cvxpy as cp


def psd_root(matrix,tolerance=1e-9):
    matrix=np.asarray(matrix,float)
    if matrix.ndim!=2 or matrix.shape[0]!=matrix.shape[1] or not np.isfinite(matrix).all():
        raise ValueError('Finite square Gram required')
    values,vectors=np.linalg.eigh((matrix+matrix.T)/2)
    if values[0]<-tolerance*max(1.,float(values[-1])):raise ValueError('Substantially indefinite Gram')
    return np.sqrt(np.maximum(values,0.))[:,None]*vectors.T


def solve_pose_cut_problem(basis,density_gram,cubic_coefficients,supports,noise_critical,density_radius,pose_scale):
    """Optimize weights restricted to basis columns and return feasible supports.

    density_gram represents [ell, -A*basis] in the continuous density space.
    supports[j] @ w <= ||F(w)|| is required for every supplied row. Its conic
    multiplier is nonnegative, and their sum is at most pose_scale. We project
    numerical multipliers onto these constraints before exporting a valid full
    weight-space pose support. Integration pads are omitted from this inner
    model and must be charged in the independent final upper audit.
    """
    wmat=np.asarray(basis,float);c=np.asarray(cubic_coefficients,float);cuts=np.asarray(supports,float)
    if wmat.ndim!=2 or c.ndim!=2 or wmat.shape[0]!=2*c.size or wmat.shape[1]<1:
        raise ValueError('Incompatible basis and cubic coefficients')
    if cuts.ndim!=2 or cuts.shape[1]!=len(wmat) or not np.isfinite(cuts).all() or not np.isfinite(wmat).all() or not np.isfinite(c).all() or (c<0).any():
        raise ValueError('Finite supports, basis and nonnegative cubic coefficients required')
    if min(noise_critical,density_radius,pose_scale)<=0:raise ValueError('Positive objective scales required')
    if np.shape(density_gram)!=(wmat.shape[1]+1,)*2:raise ValueError('Wrong density joint Gram')
    beta=cp.Variable(wmat.shape[1]);t=cp.Variable(nonneg=True)
    noise_root=psd_root(wmat.T@wmat);density_root=psd_root(density_gram)
    n,nq=c.shape;index=np.arange(len(wmat)).reshape(n,2*nq)
    real=wmat[index[:,:nq].ravel()];imag=wmat[index[:,nq:].ravel()]
    amplitude=cp.Variable(c.size)
    amplitude_cone=cp.SOC(amplitude,cp.vstack([real@beta,imag@beta]),axis=0)
    objective=noise_critical*cp.norm(noise_root@beta,2)+density_radius*cp.norm(density_root@cp.hstack([1.,beta]),2)+pose_scale*t+c.ravel()@amplitude
    constraint=t>=cuts@wmat@beta if len(cuts) else None
    problem=cp.Problem(cp.Minimize(objective),[amplitude_cone]+([] if constraint is None else [constraint]))
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-8,tol_gap_rel=1e-8,tol_feas=1e-8,max_iter=150)
    if beta.value is None or not np.isfinite(beta.value).all():raise RuntimeError(f'No finite conic candidate: {problem.status}')
    multipliers=np.maximum(np.asarray(constraint.dual_value,float),0.) if constraint is not None else np.zeros(0)
    multipliers*=min(1.,pose_scale/max(float(multipliers.sum()),np.finfo(float).tiny))
    # At an almost-zero complex weight pair, the normalized primal gradient
    # is unstable and need not approximate the optimal subgradient. The cone
    # dual supplies the appropriate interior vector. Project each pair onto
    # its known cubic dual ball to ensure global support feasibility.
    pair_support=-np.asarray(amplitude_cone.dual_value[1],float)
    pair_support*=np.minimum(1.,c.ravel()/np.maximum(np.linalg.norm(pair_support,axis=0),np.finfo(float).tiny))[None]
    cubic_support=np.concatenate([pair_support[0].reshape(n,nq),pair_support[1].reshape(n,nq)],axis=1).ravel()
    return {'weights':wmat@beta.value,'coefficients':np.asarray(beta.value),
            'pose_support':multipliers@cuts if len(cuts) else np.zeros(len(wmat)),
            'cubic_support':cubic_support,
            'support_multipliers':multipliers,'model_objective':float(problem.value),
            'model_pose_norm':float(t.value),'status':problem.status,
            'solver_iterations':int(problem.solver_stats.num_iters)}


def solve_full_weight_cut_problem(observation_projection,target_projection,cubic_coefficients,supports,noise_critical,density_radius,pose_scale):
    """Conic optimization over every observation weight and a projected density.

    The density dual lives in the explicit continuous projected basis. Its
    norm/support constraints are projected feasible before export. The caller
    must scale all dual supports for the FULL observation-space noise constraint
    and independently audit the continuous primal residual and pose bound.
    """
    f=np.asarray(observation_projection,float);b=np.asarray(target_projection,float)
    c=np.asarray(cubic_coefficients,float);cuts=np.asarray(supports,float)
    if f.ndim!=2 or b.shape!=(f.shape[1],) or c.ndim!=2 or f.shape[0]!=2*c.size or cuts.ndim!=2 or cuts.shape[1]!=len(f):
        raise ValueError('Compatible projection, cubic and cut dimensions required')
    if not all(np.isfinite(x).all() for x in [f,b,c,cuts]) or (c<0).any() or min(noise_critical,density_radius,pose_scale)<=0:
        raise ValueError('Finite arrays and positive objective scales required')
    m=len(f);n,nq=c.shape;w=cp.Variable(m);tp=cp.Variable(nonneg=True);td=cp.Variable();amplitude=cp.Variable(c.size)
    density_cone=cp.SOC(td,b-f.T@w)
    index=np.arange(m).reshape(n,2*nq)
    amplitude_cone=cp.SOC(amplitude,cp.vstack([w[index[:,:nq].ravel()],w[index[:,nq:].ravel()]]),axis=0)
    spectral=tp>=cuts@w if len(cuts) else None
    objective=noise_critical*cp.norm(w)+density_radius*td+pose_scale*tp+c.ravel()@amplitude
    constraints=[density_cone,amplitude_cone]+([] if spectral is None else [spectral])
    problem=cp.Problem(cp.Minimize(objective),constraints)
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-7,tol_gap_rel=1e-7,tol_feas=1e-7,max_iter=120)
    if w.value is None or not np.isfinite(w.value).all():raise RuntimeError(f'No full-weight candidate: {problem.status}')
    density=-np.asarray(density_cone.dual_value[1],float).ravel()
    # This dual pairs positively with b; its stationarity term is -F*density.
    density*=min(1.,density_radius/max(np.linalg.norm(density),1e-300))
    multipliers=np.maximum(np.asarray(spectral.dual_value,float),0.) if spectral is not None else np.zeros(0)
    multipliers*=min(1.,pose_scale/max(float(multipliers.sum()),1e-300))
    gp=multipliers@cuts if len(cuts) else np.zeros(m)
    pair=-np.asarray(amplitude_cone.dual_value[1],float)
    pair*=np.minimum(1.,c.ravel()/np.maximum(np.linalg.norm(pair,axis=0),1e-300))[None]
    gr=np.concatenate([pair[0].reshape(n,nq),pair[1].reshape(n,nq)],axis=1).ravel()
    defect=float(np.linalg.norm(f@density-gp-gr));dual_scale=min(1.,noise_critical/max(defect,1e-300))
    lower=max(0.,float(dual_scale*(b@density)))
    return {'weights':np.asarray(w.value),'pose_support':gp,'cubic_support':gr,'density_dual':density,
            'full_dual_lower_bound':lower,'noise_constraint_norm_before_scaling':defect,
            'dual_scale':dual_scale,'model_objective':float(problem.value),'model_pose_norm':float(tp.value),
            'status':problem.status,'solver_iterations':int(problem.solver_stats.num_iters)}
