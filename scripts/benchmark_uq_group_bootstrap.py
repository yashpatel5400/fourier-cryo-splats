#!/usr/bin/env python3
"""Conditional fixed-pose exposure-group bootstrap versus same-estimator audits.

Noise replications are independent, while bootstrap group-count draws are shared
across them to reduce computational cost; conditional MC errors refer to that
fixed resampling plan. This does not rerun pose estimation or certify real data.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.stats import norm,beta
from fourier_splats.uq_data import particle_geometry,VoxelObservationOperator,SupportedVoxelOperator,VoxelReference,local_weights
from fourier_splats.uq_physics import density_functionals
from fourier_splats.uq_subspace import orthonormalize_columns
from fourier_splats.uncertainty import bias_aware_half_width
from audit_uq_representation import dictionary
ROOT=Path(__file__).resolve().parents[1];MAPS={'10028':'2660','10049':'6487','10076':'8434'}


def binomial_record(hit):
    k=int(np.sum(hit));n=len(hit)
    return {'covered':k,'replicates':n,'coverage':k/n,'interval_95':[float(beta.ppf(.025,k,n-k+1)) if k else 0.,
        float(beta.ppf(.975,k+1,n-k)) if k<n else 1.]}


def run(args,dataset):
    rng=np.random.default_rng(args.seed);g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=args.particles,seed=609311)
    box=args.box;grid=np.arange(-box//2,box//2)/box;z,y,x=np.meshgrid(grid,grid,grid,indexing='ij')
    xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1);mask=(xyz**2).sum(-1)<=.35**2
    ck=np.load(ROOT/'results/uncertainty/development/representation'/dataset/'real_particles-spacing-2.0.npz')
    pilot=density_functionals(xyz[mask],ck['centers'],float(ck['sigma']))@ck['coefficients'];pilot/=np.linalg.norm(pilot)
    raw=SupportedVoxelOperator(VoxelObservationOperator(g['k'],g['ctf'],box),mask)
    noise=np.sqrt(np.mean(raw.forward(pilot)**2)/.1)
    op=SupportedVoxelOperator(VoxelObservationOperator(g['k'],g['ctf'],box,noise),mask)
    # Independent exact voxel columns; at this bounded size dense Grams are practical.
    phase=np.exp(-2j*np.pi*g['k']@op.xyz.T/box)*op.op.transfer[...,None]
    ambient=np.concatenate([phase.real,phase.imag],axis=1).reshape(op.shape);del phase
    ref=VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map',box=box).volume.ravel()[mask];ref/=np.linalg.norm(ref)
    delta=ref-pilot;B=2.;signal=ambient@delta
    ell=np.stack([local_weights(box,[0,0,0],.07*box).ravel()[mask],
                  (local_weights(box,[0,0,.08*box],.07*box)-local_weights(box,[0,0,-.08*box],.07*box)).ravel()[mask]],axis=1)
    observed=signal[:,None]+rng.normal(size=(len(signal),args.replicates))
    labels,groups=np.unique(g['groups'],return_inverse=True);ng=len(labels);n=op.op.n;m=2*op.op.q
    counts=rng.multinomial(ng,np.ones(ng)/ng,size=args.bootstrap)
    sizes=np.bincount(groups,minlength=ng);records=[]
    out=ROOT/'results/uncertainty/development'/args.output
    for model in args.models.split(','):
        begin=time.perf_counter()
        if model=='ambient':
            a=ambient;l=ell;s=None
        else:
            spacing={'gaussian33':3.,'gaussian123':2.}[model]
            sigma={'gaussian33':1.2,'gaussian123':1.}[model]
            s=orthonormalize_columns(density_functionals(xyz[mask],dictionary(6.,spacing),sigma))
            a=ambient@s;l=s.T@ell
        p=a.shape[1];grams=np.empty((ng,p,p));rhs=np.empty((ng,p,args.replicates))
        ab=a.reshape(n,m,p);yb=observed.reshape(n,m,args.replicates)
        for idx in range(ng):
            ag=ab[groups==idx].reshape(-1,p);yg=yb[groups==idx].reshape(-1,args.replicates)
            grams[idx]=ag.T@ag;rhs[idx]=ag.T@yg
        normal=grams.sum(axis=0);ridge=args.regularization*np.median(np.diag(normal));normal.flat[::p+1]+=ridge
        factor=cho_factor(normal,check_finite=False);v=cho_solve(factor,l,check_finite=False);w=a@v
        estimates=w.T@observed;truth=ell.T@delta;sd=np.linalg.norm(w,axis=0)
        bias=w.T@signal-truth;post_half=norm.ppf(.975)*np.sqrt(np.sum(l*v,axis=0))
        audited_bias=B*np.linalg.norm(ell-ambient.T@w,axis=0)
        audited_half=np.array([bias_aware_half_width(sv,b) for sv,b in zip(sd,audited_bias)])
        boot=np.empty((args.bootstrap,2,args.replicates))
        for idx,multiplicity in enumerate(counts):
            normal=np.tensordot(multiplicity,grams,axes=(0,0));normal.flat[::p+1]+=ridge*(multiplicity@sizes)/n
            vb=cho_solve(cho_factor(normal,check_finite=False,overwrite_a=True),l,check_finite=False)
            boot[idx]=vb.T@np.tensordot(multiplicity,rhs,axes=(0,0))
            if (idx+1)%25==0:print(dataset,model,'bootstrap',idx+1,'/',args.bootstrap,flush=True)
        lower,upper=np.quantile(boot,[.025,.975],axis=0)
        methods={'percentile':(lower,upper),'basic':(2*estimates-upper,2*estimates-lower)}
        for name,half in [('sampling_covariance',norm.ppf(.975)*sd),('Gaussian_posterior',post_half),('same_estimator_ambient_audit',audited_half)]:
            methods[name]=(estimates-half[:,None],estimates+half[:,None])
        summaries=[]
        for idx,target in enumerate(['center','contrast']):
            for name,(lo,hi) in methods.items():
                record={'target':target,'method':name,**binomial_record((lo[idx]<=truth[idx])&(hi[idx]>=truth[idx])),
                    'median_width_fraction_of_no_data':float(np.median(hi[idx]-lo[idx])/(2*B*np.linalg.norm(ell[:,idx]))),
                    'true_increment':float(truth[idx]),'estimator_bias':float(bias[idx]),'noise_sd':float(sd[idx])}
                if name not in ['percentile','basic']:
                    half=(hi[idx,0]-lo[idx,0])/2
                    record['analytic_coverage']=float(norm.cdf((half-bias[idx])/sd[idx])-norm.cdf((-half-bias[idx])/sd[idx]))
                summaries.append(record)
        records.append({'model':model,'parameters':p,'ridge':float(ridge),'group_gram_GB':float(grams.nbytes/1e9),
                        'seconds':time.perf_counter()-begin,'results':summaries})
        np.savez(out/f'{dataset}-{model}-draws.npz',bootstrap=boot,estimates=estimates,truth=truth,group_counts=counts)
        (out/f'{dataset}.json').write_text(json.dumps({'stage':'development fixed-pose known Gaussian simulation; group bootstrap does not refit poses',
            'dataset':dataset,'config':vars(args),'source_groups':ng,'particles':n,'noise_std':float(noise),
            'resampling_note':'Group-count plan reused across independent measurement-noise replicates; binomial intervals are conditional on this plan.',
            'regularization_note':'Ridge multiplier fixed before inference; bootstrap ridge scales with resampled particle count.',
            'reference_pilot_distance':float(np.linalg.norm(delta)),'models':records},indent=2)+'\n')
        print(dataset,model,'finished',records[-1]['seconds'],flush=True)
        del grams,rhs,boot,normal,factor

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--box',type=int,default=24)
    p.add_argument('--particles',type=int,default=128);p.add_argument('--replicates',type=int,default=200);p.add_argument('--bootstrap',type=int,default=199)
    p.add_argument('--models',default='gaussian33,gaussian123,ambient');p.add_argument('--regularization',type=float,default=.1)
    p.add_argument('--seed',type=int,default=609336);p.add_argument('--output',default='group-bootstrap');args=p.parse_args()
    (ROOT/'results/uncertainty/development'/args.output).mkdir(parents=True,exist_ok=True)
    for dataset in args.datasets.split(','):run(args,dataset)
