#!/usr/bin/env python3
"""Audit and enrich Gaussian-fit confidence bounds on an ambient voxel grid.

All outcomes are fixed-pose development results. The larger voxel model has a
predeclared support; it is not a guarantee over arbitrary continuous densities.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import VoxelObservationOperator,particle_geometry,local_weights,VoxelReference
from fourier_splats.uq_physics import density_functionals
from fourier_splats.uq_subspace import enrich_certificate,orthonormalize_columns,matrix_free_certificate
from fourier_splats.uncertainty import optimize_certificate,bias_aware_half_width
from audit_uq_representation import dictionary

ROOT=Path(__file__).resolve().parents[1]
MAPS={'10028':'2660','10049':'6487','10076':'8434'}


class SupportedOperator:
    def __init__(self,op,mask):self.op=op;self.mask=mask.ravel();self.shape=op.shape
    def forward(self,v):return self.op.forward(v*self.mask)
    def adjoint(self,w):return self.op.adjoint(w)*self.mask
    def forward_columns(self,s):return self.op.forward_columns(s*self.mask[:,None])


def run(args,dataset):
    t=time.perf_counter();g=particle_geometry(ROOT,dataset,'inference_half0',radius=args.frequency_radius,count=args.particles,seed=args.seed)
    box=args.box;grid=np.arange(-box//2,box//2)/box;z,y,x=np.meshgrid(grid,grid,grid,indexing='ij')
    xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    mask=(np.sum(xyz**2,axis=1)<=args.support**2).astype(float)
    pilot_path=ROOT/'results/uncertainty/development/representation'/dataset/'real_particles-spacing-2.0.npz'
    ck=np.load(pilot_path)
    pilot=density_functionals(xyz,ck['centers'],float(ck['sigma']))@ck['coefficients']*mask
    pilot/=np.linalg.norm(pilot)
    noise_free=VoxelObservationOperator(g['k'],g['ctf'],box)
    pilot_images=noise_free.forward(pilot)
    noise=np.sqrt(np.mean(pilot_images**2)/args.snr)
    op=SupportedOperator(VoxelObservationOperator(g['k'],g['ctf'],box,noise),mask)
    centers=dictionary(args.dictionary_radius,args.spacing)
    initial=density_functionals(xyz,centers,args.sigma)*mask[:,None]
    S=orthonormalize_columns(initial);a=op.forward_columns(S)
    ref=VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map',box=box)
    reference=ref.volume.ravel()*mask;reference/=np.linalg.norm(reference)
    result={'stage':'development; fixed published poses/CTFs; explicit supported finite voxel model',
            'dataset':dataset,'config':vars(args),'particles_used':len(g['indices']),'active_voxels':int(mask.sum()),
            'initial_dictionary_rank':S.shape[1],'pilot_sha256':hashlib.sha256(pilot_path.read_bytes()).hexdigest(),
            'noise_std':float(noise),'reference_pilot_distance':float(np.linalg.norm(reference-pilot)),
            'setup_seconds':time.perf_counter()-t,'targets':[]}
    for width in map(float,args.widths.split(',')):
        for label,locs,signs in [('center',[[0,0,0]],[1]),('z_contrast',[[0,0,.08],[0,0,-.08]],[1,-1])]:
            ell=sum(sign*local_weights(box,np.array(center)*box,width*box).ravel() for center,sign in zip(locs,signs))*mask
            B=args.radius;no_data=B*np.linalg.norm(ell)
            reduced=optimize_certificate(a,np.zeros((1,len(a),0)),S.T@ell,B,maxiter=300,rtol=.0005)
            full_bias=B*np.linalg.norm(ell-op.adjoint(reduced.weights))
            methods={'restricted_dictionary':(reduced.weights,reduced.half_width),
                     'audited_dictionary':(reduced.weights,bias_aware_half_width(reduced.noise_sd,full_bias))}
            begin=time.perf_counter();fit=enrich_certificate(op,ell,S,B,max_rounds=args.rounds,rtol=args.tolerance)
            methods['enriched_ambient']=(fit.weights,fit.half_width)
            enrichment_seconds=time.perf_counter()-begin
            full_record={}
            if args.matrix_free:
                begin_full=time.perf_counter();full=matrix_free_certificate(op,ell,B,rtol=args.tolerance)
                methods['matrix_free_ambient']=(full.weights,full.half_width)
                full_record={'half_width_fraction':full.half_width/no_data,'converged':full.converged,
                             'history':full.history,'seconds':time.perf_counter()-begin_full}
                print(dataset,label,width,'matrix-free',full_record['seconds'],full.history[-1],flush=True)
            rows=[]
            for name,(w,_) in methods.items():
                delta=ell-op.adjoint(w);delta*=B/max(np.linalg.norm(delta),1e-300)
                for estimator,(weights,half) in methods.items():
                    bias=weights@op.forward(delta)-ell@delta;sd=np.linalg.norm(weights)
                    coverage=float(norm.cdf((half-bias)/sd)-norm.cdf((-half-bias)/sd)) if sd>1e-15 else float(abs(bias)<=half+1e-12)
                    rows.append({'truth':'ambient_adversary_to_'+name,'method':estimator,'analytic_coverage':coverage,'bias':float(bias)})
            delta=reference-pilot
            for estimator,(weights,half) in methods.items():
                bias=weights@op.forward(delta)-ell@delta;sd=np.linalg.norm(weights)
                coverage=float(norm.cdf((half-bias)/sd)-norm.cdf((-half-bias)/sd)) if sd>1e-15 else float(abs(bias)<=half+1e-12)
                rows.append({'truth':'independent_EMDB_voxel_reference','method':estimator,'analytic_coverage':coverage,
                             'bias':float(bias),'density_bound_satisfied':bool(np.linalg.norm(delta)<=B)})
            record={'target':label,'width_fraction_field':width,'width_A':width*g['field_A'],
                    'no_data_half_width':float(no_data),'restricted_width_fraction':float(reduced.half_width/no_data),
                    'audited_width_fraction':float(methods['audited_dictionary'][1]/no_data),
                    'enriched_width_fraction':float(fit.half_width/no_data),'enriched_dimension':fit.basis_dimension,
                    'converged':fit.converged,'seconds':enrichment_seconds,'history':fit.history,'coverage':rows,
                    'matrix_free':full_record}
            result['targets'].append(record)
            print(dataset,label,width,{k:v for k,v in record.items() if k not in ['history','coverage','matrix_free']},flush=True)
    result['seconds']=time.perf_counter()-t
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076')
    p.add_argument('--particles',type=int,default=128);p.add_argument('--box',type=int,default=24)
    p.add_argument('--frequency-radius',type=float,default=5.);p.add_argument('--support',type=float,default=.35)
    p.add_argument('--dictionary-radius',type=float,default=6.);p.add_argument('--spacing',type=float,default=3.)
    p.add_argument('--sigma',type=float,default=1.2);p.add_argument('--widths',default='.03,.07')
    p.add_argument('--snr',type=float,default=.1);p.add_argument('--radius',type=float,default=2.)
    p.add_argument('--rounds',type=int,default=40);p.add_argument('--tolerance',type=float,default=.005)
    p.add_argument('--matrix-free',action='store_true')
    p.add_argument('--seed',type=int,default=609311);p.add_argument('--output',default='ambient-audit.json');args=p.parse_args()
    output=ROOT/'results/uncertainty/development'/args.output;result={'config':vars(args),'cases':[]}
    for dataset in args.datasets.split(','):
        result['cases'].append(run(args,dataset));output.write_text(json.dumps(result,indent=2)+'\n')
