#!/usr/bin/env python3
"""Postprocess a completed estimator's dual certificate without refitting it."""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh, ArpackNoConvergence
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_pose_optimization import cubic_penalty_coefficients, amplitude_gradient
from fourier_splats.uq_dual_mixture import closest_support_mixture
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--threads',type=int,default=1)
    p.add_argument('--modes',type=int,default=12);args=p.parse_args()
    if args.threads>1 and sys.platform=='darwin':
        import finufft
        if 'torch' in sys.modules or 'tmp/nufft-openmp' not in finufft.__file__:
            raise RuntimeError('Use isolated CPU-only OpenMP runtime')
    source=ROOT/args.fit;prior=json.loads(source.read_text())
    if not prior.get('complete') or 'fit' not in prior:raise RuntimeError('Completed fit required')
    dataset=prior['dataset'];target=prior['target'];degrees=prior['rotation_radius_degrees']
    out=BASE/'pose-dual-mixture'/source.parent.name;out.mkdir(parents=True,exist_ok=True);path=out/source.name
    if path.exists():raise RuntimeError('Preserve previous dual audit')
    start=time.perf_counter();snapshot=source_snapshot(ROOT,Path(__file__))
    result={'stage':'Exploratory dual-only postprocessing; the fitted estimator and its primal interval are unchanged',
            'config':vars(args),'source_snapshot':snapshot,'source_fit_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'dataset':dataset,'target':target,'rotation_radius_degrees':degrees,'complete':False}
    def save():path.write_text(json.dumps(result,indent=2)+'\n')
    save()
    try:
        cfg=prior['source_geometry_config'];g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=128,seed=cfg['seed'])
        saved=np.load(source.with_name(source.stem+'-weights.npz'));w=saved['weights'];noise=float(saved['noise_std'])
        initial=np.load(BASE/'continuous-quadrature-optimized'/f'{dataset}-{target}-0.07-weights.npz')['weights']
        np.testing.assert_array_equal(saved['indices'],g['indices'])
        angle=np.deg2rad(degrees);shift=prior['translation_radius_A']/g['field_A'];fit=prior['fit']
        op=PolynomialPoseFieldOperator(g['k'],g['q'],g['ctf'],initial,noise,angle,shift,order=fit['pose_quadrature_order'],nthreads=args.threads)
        scaling=op.establish_group_scaling();scale_sum=float(np.where(scaling['group_scales']>0,scaling['group_scales'],1.).sum())
        np.testing.assert_allclose(scale_sum,fit['initial_group_scale_sum'],rtol=1e-8)
        op.set_weights(w);operator=LinearOperator((op.shape[0],)*2,matvec=op.spatial_gram,dtype=float)
        mode_status='all requested modes converged'
        try:
            values,vectors=eigsh(operator,k=args.modes,which='LA',tol=1e-4,maxiter=100,
                                 ncv=max(25,2*args.modes+1),v0=np.random.default_rng(609911).normal(size=op.shape[0]))
        except ArpackNoConvergence as exc:
            if exc.eigenvalues is None or len(exc.eigenvalues)==0:raise
            values,vectors=exc.eigenvalues,exc.eigenvectors;mode_status=str(exc)
        vectors/=np.linalg.norm(vectors,axis=0)
        directions=[]
        for v in vectors.T:
            u=op.rmatvec(v);directions.append(3*np.sqrt(scale_sum)*op.weight_gradient(u/max(np.linalg.norm(u),1e-300),v))
        directions=np.asarray(directions)
        gram=QuadratureObservationGram(g['k'],g['ctf'],noise,order=40,preconditioner_rank=0);gram.nthreads=args.threads
        centers=[[0,0,0]] if target=='center' else [[0,0,.08],[0,0,-.08]];signs=[1] if target=='center' else [1,-1]
        a,norm2=gram.target(centers,signs,.07);gw=gram.matvec(w);error=gram.quadrature_error(w)
        residual=np.sqrt(max(0,norm2-2*w@a+w@gw+error['squared_field_norm']))
        factor=2/max(residual,1e-300)
        remainder=cubic_penalty_coefficients(g['k'],g['q'],g['ctf']/noise,angle,shift,3.)
        gr=amplitude_gradient(remainder,w)
        mixture=closest_support_mixture(directions,factor*(a-gw)-gr)
        beta=mixture.pop('mixture');mixture['mixture']=beta.tolist()
        defect=mixture['distance']+factor*error['gram_action_norm']
        z=norm.isf(fit['alpha_noise']/2);dual_scale=min(1.,z/max(defect,1e-300))
        candidate=max(0.,float(dual_scale*factor*(norm2-w@a)))
        lower=max(fit['dual_lower_bound'],candidate);upper=fit['sum_objective_upper']
        result.update(complete=True,seconds=time.perf_counter()-start,available_modes=len(values),eigensolver_status=mode_status,
            mixture_fit=mixture,candidate_dual_lower_bound=candidate,retained_dual_lower_bound=lower,
            original_dual_lower_bound=fit['dual_lower_bound'],primal_upper_unchanged=upper,
            original_relative_gap=fit['relative_sum_gap'],tightened_relative_gap=float((upper-lower)/upper),
            numerical_failure=bool(lower>upper*(1+1e-6)))
        save();print('DONE',result['original_relative_gap'],result['tightened_relative_gap'],flush=True)
        if result['numerical_failure']:raise AssertionError('Dual exceeds primal; full diagnostic retained')
    except Exception as exc:
        result.update(error=repr(exc),seconds=time.perf_counter()-start);save();raise


if __name__=='__main__':main()
