#!/usr/bin/env python3
"""Audit coarse-grid weights on a common finer physical-density grid.

All centers/truth evaluations use the same fine grid. Discretization is not
silently changed when reporting coverage. Unit L2 coordinates include the
voxel-volume factor on the common unit field.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry,VoxelObservationOperator,SupportedVoxelOperator,local_weights,VoxelReference
from fourier_splats.uq_physics import density_functionals
from fourier_splats.uq_subspace import matrix_free_certificate
from fourier_splats.uncertainty import bias_aware_half_width
ROOT=Path(__file__).resolve().parents[1]
MAPS={'10028':'2660','10049':'6487','10076':'8434'}


def model(g,ck,box,noise=None):
    grid=np.arange(-box//2,box//2)/box;z,y,x=np.meshgrid(grid,grid,grid,indexing='ij')
    xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1);mask=(xyz**2).sum(axis=-1)<=.35**2
    pilot=density_functionals(xyz[mask],ck['centers'],float(ck['sigma']))@ck['coefficients'];pilot/=np.linalg.norm(pilot)
    if noise is None:
        raw=SupportedVoxelOperator(VoxelObservationOperator(g['k'],g['ctf'],box,box**1.5),mask)
        noise=np.sqrt(np.mean(raw.forward(pilot)**2)/.1)
    op=SupportedVoxelOperator(VoxelObservationOperator(g['k'],g['ctf'],box,noise*box**1.5),mask)
    return op,pilot,mask,noise


def target(box,mask,width,name):
    locations=[([0,0,0],1)] if name=='center' else [([0,0,.08],1),([0,0,-.08],-1)]
    return box**1.5*sum(s*local_weights(box,np.array(c)*box,width*box).ravel()[mask] for c,s in locations)


def run(args,dataset):
    g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=128,seed=609311)
    ck=np.load(ROOT/'results/uncertainty/development/representation'/dataset/'real_particles-spacing-2.0.npz')
    fine,pilot,mask,noise=model(g,ck,args.fine_box)
    ref=VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map',box=args.fine_box).volume.ravel()[mask]
    ref/=np.linalg.norm(ref);delta=ref-pilot;signal=fine.forward(delta)
    result={'dataset':dataset,'fine_box':args.fine_box,'noise_std_physical_units':float(noise),
            'reference_pilot_distance':float(np.linalg.norm(delta)),'targets':[]}
    for box in map(int,args.boxes.split(',')):
        op,_,coarse_mask,_=model(g,ck,box,noise)
        for width in [.03,.07]:
            for name in ['center','contrast']:
                ell=target(box,coarse_mask,width,name);ellfine=target(args.fine_box,mask,width,name);B=2.
                begin=time.perf_counter();fit=matrix_free_certificate(op,ell,B,rtol=.005)
                sd=np.linalg.norm(fit.weights);fullbias=B*np.linalg.norm(ellfine-fine.adjoint(fit.weights))
                qfine=bias_aware_half_width(sd,fullbias);actualbias=float(fit.weights@signal-ellfine@delta)
                analytic=lambda q:float(norm.cdf((q-actualbias)/sd)-norm.cdf((-q-actualbias)/sd)) if sd>1e-15 else float(abs(actualbias)<=q)
                record={'box':box,'width_fraction':width,'target':name,'coarse_gap_converged':fit.converged,
                        'coarse_width':fit.half_width,'fine_audited_width':qfine,'fine_no_data_width':float(B*np.linalg.norm(ellfine)),
                        'relative_audit_inflation':float(qfine/fit.half_width),
                        'coarse_reference_coverage':analytic(fit.half_width),'fine_audited_reference_coverage':analytic(qfine),
                        'fine_reference_target':float(ellfine@ref),'actual_bias':actualbias,'noise_sd':float(sd),
                        'seconds':time.perf_counter()-begin,'iterations':len(fit.history)}
                result['targets'].append(record)
                print(dataset,box,width,name,'width audit factor',record['relative_audit_inflation'],'seconds',record['seconds'],flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--boxes',default='16,24,32,48');p.add_argument('--fine-box',type=int,default=64)
    args=p.parse_args();out=ROOT/'results/uncertainty/development/grid-refinement.json'
    result={'stage':'development fixed-pose finite-grid sensitivity; all inference centers/true targets use common fine grid',
            'config':vars(args),'cases':[]}
    for dataset in args.datasets.split(','):
        result['cases'].append(run(args,dataset));out.write_text(json.dumps(result,indent=2)+'\n')
