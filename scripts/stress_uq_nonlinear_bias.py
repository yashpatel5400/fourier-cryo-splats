#!/usr/bin/env python3
"""Search the joint density/pose class for a large bias of fixed saved weights.

The density optimization is exact; pose optimization is local, so found biases
are lower bounds only. All reported candidates are recomputed in NumPy float64.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
import torch
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_nonlinear_stress import TorchPoseAdjoint,exact_pose_adjoint,eliminated_bias
from audit_uq_grid_refinement import model
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def run(args,dataset):
    folder=BASE/args.source;source=json.loads((folder/f'{dataset}.json').read_text());cfg=source['config']
    saved=np.load(folder/f'{dataset}-coarse-weights.npz');w=saved['weights'];ell=saved['ell'];pilot=saved['pilot']
    g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=cfg['particles'],seed=609315)
    np.testing.assert_array_equal(g['indices'],saved['indices'])
    ck=np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
    op,recreated,_,_=model(g,ck,cfg['coarse_box'],source['noise_std_physical'])
    np.testing.assert_allclose(pilot,recreated,atol=1e-12)
    B=2.;angle=np.deg2rad(cfg['angle']);shift=.01;n=op.op.n
    grid=next(r for r in source['grids'] if r['box']==cfg['coarse_box'])
    upper=grid['density_bias']+sum(grid['bias_components'].values());half=source['coarse_half_width'];sd=np.linalg.norm(w)
    device='mps' if args.device=='auto' and torch.backends.mps.is_available() else args.device
    if device=='auto':device='cpu'
    dtype=torch.float32 if device=='mps' else torch.float64;torch.set_num_threads(4)
    forward=TorchPoseAdjoint(op,g['q'],w,angle,shift,device,dtype)
    tensor=lambda x:torch.as_tensor(np.asarray(x).copy(),dtype=dtype,device=device)
    nominal=forward(torch.zeros((n,5),dtype=dtype,device=device)).detach();ellt=tensor(ell);pilott=tensor(pilot)
    exact_nominal=exact_pose_adjoint(op,g['q'],w,np.zeros((n,5)),angle,shift)
    rng=np.random.default_rng(args.seed);rows=[];out=BASE/args.output;out.mkdir(parents=True,exist_ok=True)
    result={'stage':'development feasible lower bounds on worst bias; local pose optimization is not a global solution',
            'dataset':dataset,'config':vars(args),'source_config':cfg,'source_json_sha256':hashlib.sha256((folder/f'{dataset}.json').read_bytes()).hexdigest(),
            'weights_sha256':source['weights_sha256'],'torch_version':torch.__version__,'device':device,'dtype':str(dtype),
            'bias_upper_certificate':float(upper),'half_width':float(half),'noise_sd':float(sd),'candidates':rows,
            'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'src/fourier_splats/uq_nonlinear_stress.py']}}
    for sign in [-1,1]:
        for restart in range(args.restarts):
            start=time.perf_counter();initial=rng.normal(size=(n,5));initial/=np.linalg.norm(initial,axis=1,keepdims=True)
            if restart==0:initial*=0
            u=tensor(initial);u.requires_grad_(True);optimizer=torch.optim.Adam([u],lr=args.learning_rate)
            best=-np.inf;best_u=None;trace=[]
            for step in range(args.steps+1):
                optimizer.zero_grad();value=eliminated_bias(forward(u),nominal,ellt,pilott,B,sign)
                scalar=float(value.detach().cpu())
                if scalar>best:best=scalar;best_u=u.detach().cpu().numpy().copy()
                if step%10==0 or step==args.steps:trace.append({'step':step,'best_bias_float':best})
                if step==args.steps:break
                (-value).backward();optimizer.step()
                with torch.no_grad():u.div_(torch.clamp(torch.linalg.vector_norm(u,dim=1,keepdim=True),min=1))
            # Float32 projection may lie just outside the ball: reproject in float64.
            best_u=best_u.astype(np.float64);best_u/=np.maximum(1,np.linalg.norm(best_u,axis=1,keepdims=True))
            exact=exact_pose_adjoint(op,g['q'],w,best_u,angle,shift);h=exact-ell
            delta=sign*B*h/np.linalg.norm(h);bias=float(pilot@(exact-exact_nominal)+delta@h)
            signed=sign*bias;coverage=float(norm.cdf((half-bias)/sd)-norm.cdf((-half-bias)/sd))
            record={'sign':sign,'restart':restart,'signed_bias_lower_bound':signed,'absolute_bias':abs(bias),
                    'fraction_of_bias_certificate':float(abs(bias)/upper),'analytic_coverage':coverage,
                    'maximum_pose_norm':float(np.linalg.norm(best_u,axis=1).max()),'density_delta_norm':float(np.linalg.norm(delta)),
                    'float_optimization_bias':best,'float_vs_double_bias_error':abs(best-signed),'seconds':time.perf_counter()-start,'trace':trace}
            rows.append(record);np.savez(out/f'{dataset}-sign{sign}-restart{restart}.npz',pose=best_u,density_delta=delta)
            result['maximum_found_bias_fraction']=max(r['fraction_of_bias_certificate'] for r in rows)
            result['minimum_candidate_coverage']=min(r['analytic_coverage'] for r in rows)
            (out/f'{dataset}.json').write_text(json.dumps(result,indent=2)+'\n')
            print(dataset,sign,restart,'bias/upper',record['fraction_of_bias_certificate'],'coverage',coverage,'seconds',record['seconds'],flush=True)
            if abs(bias)>upper+1e-7*max(1,upper):raise AssertionError('Feasible adversary exceeds certified upper bound')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',default='pose-grid-center-0.5');p.add_argument('--output',default='nonlinear-bias-stress')
    p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--steps',type=int,default=100);p.add_argument('--restarts',type=int,default=3)
    p.add_argument('--device',choices=['cpu','mps','auto'],default='auto');p.add_argument('--learning-rate',type=float,default=.05);p.add_argument('--seed',type=int,default=609352)
    args=p.parse_args()
    for dataset in args.datasets.split(','):run(args,dataset)
