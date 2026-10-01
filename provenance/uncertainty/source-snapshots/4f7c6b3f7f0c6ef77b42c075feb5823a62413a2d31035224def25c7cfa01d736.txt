"""Joint density/pose norm design using classical robust least-squares cuts.

Quadrature witnesses and numerical conic solves guide weight selection only.
A fresh, separately randomized continuous upper audit is required afterward.
No routine here claims full-space optimality or calibrated experimental inputs.
"""
import time
from numbers import Integral
from scipy.stats import norm
import numpy as np
from .uq_continuous_quadrature import QuadratureObservationGram
from .uq_cubic_subspace import ReducedCubicDesign, psd_root
from .uq_trust_region import dense_quadratic_witness, krylov_cache, quadratic_dual_upper
from .uq_cubic_pose import cubic_residual_cross
from .uq_joint_bias import joint_density_pose_bias
from .uq_random_spectral import gaussian_power_upper
from .uq_intervals import bias_aware_half_width_stable


class ReducedJointCubicDesign:
    def __init__(self, objective, basis, centers, signs, width):
        self.obj = objective; self.reduced = ReducedCubicDesign(objective,basis)
        self.Q = self.reduced.Q; self.p = self.Q.shape[1]
        op = objective.op
        gram = QuadratureObservationGram(op.k,op.ctf,objective.gram.noise_std,order=op.order,preconditioner_rank=0)
        gram.nthreads = op.nthreads
        self.nominal = np.column_stack([op.sqrt_quad*gram.field(q) for q in self.Q.T])
        target = np.zeros(len(op.xyz))
        for mu,sign in zip(np.asarray(centers,float),np.asarray(signs,float)):
            target += sign*np.exp(-np.sum((op.xyz-mu)**2,axis=1)/(2*width**2))/(2*np.pi*width**2)**1.5
        self.target = op.sqrt_quad*target
        self.roots = []; self.witnesses = []; self.cut_diagnostics = []
        self.add_cut(np.zeros(op.shape[1]))

    def add_cut(self, witness):
        v = np.asarray(witness,float); op = self.obj.op
        if v.shape != (op.shape[1],) or not np.isfinite(v).all() or np.linalg.norm(v)>self.obj.L*(1+1e-12):
            raise ValueError('Compatible feasible pose witness required')
        if self.witnesses and min(np.linalg.norm(v-old) for old in self.witnesses)<1e-10*max(1.,self.obj.L): return False
        fields = self.nominal.copy()
        if np.linalg.norm(v):
            for j,q in enumerate(self.Q.T):
                op.set_weights(q); fields[:,j] += op.matvec(v)
        augmented = np.column_stack([fields,-self.target])
        root,record = psd_root(augmented.T@augmented)
        self.roots.append(root); self.witnesses.append(v.copy()); self.cut_diagnostics.append(record)
        return True

    def ancillary_value(self,x):
        red = self.reduced; obj = self.obj
        terms = {'noise':obj.z*np.linalg.norm(red.noise_root@x),
            'pilot':obj.L*np.linalg.norm(red.pilot_root@x),
            'remainder':sum(np.linalg.norm(root@x) for root in red.remainder_roots)+red.residual_coefficients@abs(x)}
        return {k:float(v) for k,v in terms.items()}

    def witness(self,x,*,start,steps=24,dense=False):
        op=self.obj.op;op.set_weights(self.Q@x)
        h=self.target-self.nominal@x;cross=op.rmatvec(h)
        def action(v):return op.rmatvec(op.matvec(v))
        if dense:
            matrix=np.column_stack([op.matvec(e) for e in np.eye(op.shape[1])])
            gram=matrix.T@matrix; basis=np.eye(op.shape[1]); cache={'rank':op.shape[1],'dense':True}
        else:
            basis,images,cache=krylov_cache(action,cross,steps=steps,extras=[start])
            gram=.5*(basis.T@images+images.T@basis)
        if basis.shape[1]:
            small,record=dense_quadratic_witness(gram,basis.T@cross,self.obj.L)
            v=basis@small; length=np.linalg.norm(v)
            if length: v *= min(1.,np.nextafter(self.obj.L,0.)/length)
        else:
            v=np.zeros(op.shape[1]);record={'value':0.,'hard_case':False}
        norm=float(np.linalg.norm(h-op.matvec(v)))
        witness_norm=norm
        known=[float(np.linalg.norm(root@np.r_[x,1.])) for root in self.roots]
        raw=[float(np.sqrt(max(0.,value**2-diag['diagonal_pad']*(x@x+1))))
             for value,diag in zip(known,self.cut_diagnostics)]
        if max(known)>norm:
            index=int(np.argmax(known));v=self.witnesses[index].copy();norm=known[index]
        record.update(cache=cache,joint_norm_guide=norm,existing_cut_norms=known,existing_cut_unpadded_norms=raw,
            oracle_unpadded_witness_norm=witness_norm,joint_norm_unpadded_guide=max(witness_norm,max(raw)))
        return v,record

    def solve_master(self,solver='CLARABEL'):
        import cvxpy as cp
        x=cp.Variable(self.p);absolute=cp.Variable(self.p,nonneg=True);joint=cp.Variable(nonneg=True)
        constraints=[absolute>=x,absolute>=-x]
        for root in self.roots:constraints.append(cp.norm(root@cp.hstack([x,1.]))<=joint)
        red=self.reduced;obj=self.obj
        cost=obj.z*cp.norm(red.noise_root@x)+obj.B*joint+obj.L*cp.norm(red.pilot_root@x)
        cost+=sum(cp.norm(root@x) for root in red.remainder_roots)+red.residual_coefficients@absolute
        problem=cp.Problem(cp.Minimize(cost),constraints);begin=time.perf_counter()
        opts=dict(tol_gap_abs=1e-8,tol_gap_rel=1e-8,tol_feas=1e-8,max_iter=500) if solver=='CLARABEL' else dict(eps=1e-7,max_iters=100000)
        value=problem.solve(solver=solver,**opts)
        if problem.status not in ('optimal','optimal_inaccurate') or x.value is None or not np.isfinite(x.value).all():
            raise RuntimeError(f'Joint restricted master failed: {problem.status}')
        return x.value,{'solver':solver,'status':problem.status,'objective':float(value),
            'coefficients':x.value.tolist(),'seconds':time.perf_counter()-begin,'iterations':problem.solver_stats.num_iters,
            'scope':'Numerical restricted quadrature problem; not a continuous full-space lower certificate.'}


def validate_joint_configuration(objective,optimization_seed,certificate_seed,alpha,delta):
    if any(isinstance(seed,bool) or not isinstance(seed,Integral) or seed<0 for seed in [optimization_seed,certificate_seed]):
        raise ValueError('Explicit nonnegative integer seeds required')
    if optimization_seed==certificate_seed:raise ValueError('Distinct design and certificate seeds required')
    if not 0<delta<alpha<.5 or abs(objective.z-norm.isf((alpha-delta)/2))>1e-12:
        raise ValueError('Inconsistent critical value or error budget')


def rescore_joint_candidates(history,roots,density_radius,master_objective):
    """Revisit every evaluated coefficient vector against the final cut set."""
    rows=[]
    for record in history:
        x=np.asarray(record['coefficients'])
        joint=max(record['oracle']['joint_norm_guide'],max(np.linalg.norm(root@np.r_[x,1.]) for root in roots))
        ancillary=sum(value for key,value in record['terms'].items() if key!='joint_density_pose')
        value=float(ancillary+density_radius*joint)
        rows.append({'evaluation':record['evaluation'],'final_cut_guide':value,'joint_norm_guide':float(joint),
            'restricted_guide_gap':max(0.,value-master_objective)/max(value,np.finfo(float).tiny)})
    return rows,min(rows,key=lambda row:row['final_cut_guide'])


def design_joint_cubic(objective,basis,initial_weights,centers,signs,width,*,optimization_seed,certificate_seed,alpha=.05/12,delta=1e-6/12,
                       evaluations=12,krylov_steps=24,design_seconds=5400,separation_tolerance=.001,
                       dense_oracle=False,callback=None,checkpoint_callback=None):
    validate_joint_configuration(objective,optimization_seed,certificate_seed,alpha,delta)
    if evaluations<1 or krylov_steps<1 or design_seconds<=0 or not 0<=separation_tolerance<1:
        raise ValueError('Positive budgets and valid separation tolerance required')
    begin=time.perf_counter();red=ReducedJointCubicDesign(objective,basis,centers,signs,width)
    x=red.Q.T@initial_weights
    np.testing.assert_allclose(red.Q@x,initial_weights,rtol=1e-9,atol=1e-10)
    start=np.random.default_rng(optimization_seed).normal(size=objective.op.shape[1]);start/=np.linalg.norm(start)
    best={'value':np.inf};history=[];master=None;checkpointed=set()
    for index in range(evaluations):
        if history and time.perf_counter()-begin>=design_seconds:break
        stamp=time.perf_counter();v,oracle=red.witness(x,start=start,steps=krylov_steps,dense=dense_oracle)
        terms=red.ancillary_value(x);terms['joint_density_pose']=objective.B*oracle['joint_norm_guide']
        value=float(sum(terms.values()));record={'evaluation':index+1,'guide_objective':value,'terms':terms,
            'oracle':oracle,'coefficients':x.tolist(),'previous_master':master,'elapsed_seconds':time.perf_counter()-begin}
        if value<best['value']:
            best={'value':value,'weights':red.Q@x,'evaluation':index+1}
            if checkpoint_callback:checkpoint_callback(best['weights'],record);checkpointed.add(index+1)
        red.add_cut(v)
        proposal,master=red.solve_master();record['next_master']=master
        gap=max(0.,value-master['objective'])/max(value,np.finfo(float).tiny)
        record.update(restricted_guide_gap=gap,seconds=time.perf_counter()-stamp,cut_count=len(red.roots))
        history.append(record)
        if callback:callback(record)
        if index>=1 and gap<=separation_tolerance:break
        x=proposal
    if not history:raise RuntimeError('No joint candidate evaluated')
    rescored,selected=rescore_joint_candidates(history,red.roots,objective.B,master['objective'])
    selected_record=next(row for row in history if row['evaluation']==selected['evaluation'])
    selected_weights=red.Q@np.asarray(selected_record['coefficients'])
    if checkpoint_callback and selected['evaluation'] not in checkpointed:
        checkpoint_callback(selected_weights,dict(selected_record,selected_by_final_rescore=True))
    return {'weights':selected_weights,'selected_evaluation':selected['evaluation'],'selected_guide':selected['final_cut_guide'],
        'history':history,'final_cut_rescoring':rescored,'optimization_seconds':time.perf_counter()-begin,
        'optimization_seed':optimization_seed,'certificate_seed':certificate_seed,
        'restricted_guide_gap':selected['restricted_guide_gap'],'last_iteration_guide_gap':history[-1]['restricted_guide_gap'],
        'reduced_diagnostics':red.reduced.diagnostics,'cut_diagnostics':red.cut_diagnostics,
        'scope':'Final-cut-rescored joint quadrature guide only; no continuous/full-space convergence certificate.'}


def joint_upper_components(objective,weights,spectral_upper):
    """Reuse frozen upper primitives without computing an unrelated dual gap."""
    objective.op.set_weights(weights)
    residual,_,guard=objective.density(weights,final_guard=True)
    pad,_=objective.mass_terms(weights);field=objective.L*np.sqrt(spectral_upper+pad)
    pilot=objective.L*np.linalg.norm(objective.op.pair_moments(objective.moments))
    rem,_=objective.remainder.value_gradient(weights)
    bias=objective.B*(residual+field)+pilot+(objective.B+objective.P)*rem
    sd=float(np.linalg.norm(weights))
    return {'noise_sd':sd,'bias':float(bias),'density_bias':float(objective.B*residual),
        'pose_field_upper':float(field),'pose_polynomial_bias':float(objective.B*field),
        'pilot_polynomial_bias':float(pilot),'quartic_bias':float((objective.B+objective.P)*rem),
        'pose_integration_pad':pad,'density_majorant_eta':objective.eta,'density_squared_roundoff_guard':float(guard),
        'target_norm':float(np.sqrt(objective.target_norm2)),'sum_objective_upper':float(objective.z*sd+bias)}


def audit_joint_cubic(objective,weights,centers,signs,width,*,certificate_seed,optimization_seed,delta=1e-6/12,
                      alpha=.05/12,probes=4,power_iterations=40,krylov_steps=40,callback=None):
    """Continuous joint upper bound on one fresh fixed-weight spectral event."""
    validate_joint_configuration(objective,optimization_seed,certificate_seed,alpha,delta)
    begin=time.perf_counter();op=objective.op;op.set_weights(weights)
    spectral=gaussian_power_upper(op.spatial_gram,op.shape[0],delta,probes,power_iterations,certificate_seed,callback=callback)
    # Reuse the existing density/pose/remainder enclosures, but DO NOT report
    # their triangle-objective dual gap as a gap for the different joint design.
    old=joint_upper_components(objective,weights,spectral['eigenvalue_upper'])
    cross=cubic_residual_cross(op,weights,objective.gram.noise_std,centers,signs,width)
    b=cross.pop('cross_vector')
    q,gq,cache=krylov_cache(lambda v:op.rmatvec(op.matvec(v)),b,steps=krylov_steps)
    dual=quadratic_dual_upper(b,objective.L,spectral['eigenvalue_upper'],q,gq)
    h=old['density_bias']/objective.B;L=objective.L
    squared=h*h+dual['quadratic_upper']+L*L*old['pose_integration_pad']+2*L*cross['quadrature_norm_pad']
    prior=joint_density_pose_bias(h,old['pose_field_upper'],cross['norm_upper'],L,objective.B,objective.P,
        old['quartic_bias']/(objective.B+objective.P),float(np.linalg.norm(op.pair_moments(objective.moments))))
    resolved=objective.B*np.sqrt(max(0.,squared))+prior['pilot_polynomial_bias_upper']+old['quartic_bias']
    bias=min(float(resolved),prior['bias_upper'])
    return {'bias_upper':bias,'resolvent_bias_upper':float(resolved),'old_joint':prior,
        'old_triangle_audit':old,'spectral_upper_certificate':spectral,'cross':cross,
        'residual_cross_vector':b,'resolvent':dual,'krylov':cache,'joint_squared_upper':float(squared),
        'noise_sd':old['noise_sd'],'target_norm':old['target_norm'],
        'half_width':bias_aware_half_width_stable(old['noise_sd'],bias,alpha-delta),
        'alpha_total':alpha,'alpha_noise':alpha-delta,'alpha_numerical':delta,
        'certificate_seconds':time.perf_counter()-begin,'certificate_seed':certificate_seed,
        'scope':'One supplied noise/class/pose model; common fresh spectral event; real-arithmetic bounds with heuristic numerical guards.'}
