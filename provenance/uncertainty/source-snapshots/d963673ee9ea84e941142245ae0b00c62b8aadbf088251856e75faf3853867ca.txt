#!/usr/bin/env python3
"""Independent small conic optimum checks for the new pose-aware objective."""
import json
import time
from pathlib import Path
import numpy as np
import cvxpy as cp
from scipy.stats import norm
from fourier_splats.uq_continuous import continuous_certificate,ContinuousObservationGram
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_pose_optimization import pose_aware_certificate,cubic_penalty_coefficients
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/uncertainty/development/pose-optimizer-conic.json'
    if out.exists():raise RuntimeError('Preserve prior independent validation')
    snapshot=source_snapshot(ROOT,Path(__file__));start=time.perf_counter()
    q=np.array([[[.8,.2],[.1,.7]]]);k=np.pad(q,((0,0),(0,0),(0,1)))
    ctf=np.array([[.8,-.6]]);noise=.15;B=1.;P=0.;angle=.04;shift=.002;delta=1e-6
    centers=[[0,0,0]];signs=[1];width=.25
    gram=QuadratureObservationGram(k,ctf,noise,order=16,preconditioner_rank=4)
    init=continuous_certificate(gram,centers,signs,width,B,rtol=1e-5,maxiter=200)
    w0=init['weights'];m=len(w0)
    op=PolynomialPoseFieldOperator(k,q,ctf,w0,noise,angle,shift,order=8,backend='direct')
    scaling=op.establish_group_scaling();coeff=op.coefficient_bound_matrix()[0]
    fields=[]
    for e in np.eye(m):
        op.set_weights(e)
        fields.append(np.column_stack([op.matvec(u) for u in np.eye(op.shape[1])]))
    all_fields=np.concatenate(fields,axis=1);basis,compressed=np.linalg.qr(all_fields,mode='reduced')
    projection_error=float(np.linalg.norm(basis@compressed-all_fields)/np.linalg.norm(all_fields))
    a,tnorm2=gram.target(centers,signs,width)
    exact=ContinuousObservationGram(k,ctf,noise)
    g=np.column_stack([exact.matvec(e) for e in np.eye(m)])
    block=np.block([[np.array([[tnorm2]]),-a[None]],[-a[:,None],g]])
    eigen,vectors=np.linalg.eigh(block)
    if eigen.min()<-1e-10:raise AssertionError('Invalid exact density Gram')
    root=np.sqrt(np.maximum(eigen,0))[:,None]*vectors.T
    w=cp.Variable(m);amps=cp.norm(cp.vstack([w[:2],w[2:]]),axis=0)
    field=sum(w[j]*compressed[:,j*20:(j+1)*20] for j in range(m))
    spectral=cp.norm(field,2)
    pad=np.sqrt(scaling['kernel_error'])*cp.norm(coeff@amps,2)
    pose=(B+P)*np.sqrt(scaling['sum_group_scales'])*cp.norm(cp.hstack([spectral,pad]),2)
    remainder_coeff=cubic_penalty_coefficients(k,q,ctf/noise,angle,shift,B+P)[0]
    z=norm.isf((.05-delta)/2)
    objective=z*cp.norm(w,2)+B*cp.norm(root@cp.hstack([1,w]),2)+pose+remainder_coeff@amps
    problem=cp.Problem(cp.Minimize(objective))
    if not problem.is_dcp():raise AssertionError('Conic formulation must be convex')
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-8,tol_feas=1e-8,tol_gap_rel=1e-8,max_iter=200)
    if problem.status!='optimal':raise AssertionError('Independent conic optimum unavailable')
    fit=pose_aware_certificate(gram,q,centers,signs,width,B,P,w0,angle,shift,maxiter=150,
                               ritz_tolerance=1e-9,spectral_smoothing=1e-5,power_iterations=100,
                               spectral_modes=8,pose_order=8,
                               callback=lambda r:print(r,flush=True))
    fit.pop('weights')
    result={'stage':'Independent small conic verification of shared-field objective and dual lower bound',
            'source_snapshot':snapshot,'conic_solver':'CLARABEL','conic_status':problem.status,
            'conic_optimum':float(problem.value),'row_compression_relative_error':projection_error,
            'pose_order':8,'fit':fit,'seconds':time.perf_counter()-start,
            'dual_below_conic':bool(fit['dual_lower_bound']<=problem.value+1e-6),
            'conic_below_reported_upper':bool(problem.value<=fit['sum_objective_upper']+1e-6),
            'upper_over_conic':float(fit['sum_objective_upper']/problem.value)}
    out.write_text(json.dumps(result,indent=2)+'\n')
    if not result['dual_below_conic'] or not result['conic_below_reported_upper']:
        raise AssertionError('Primal/dual bracketing failed; full failure record retained')
    print('DONE',result['upper_over_conic'],fit['relative_sum_gap'],flush=True)


if __name__=='__main__':main()
