#!/usr/bin/env python3
"""Audit coarse optimized nonlinear certificates on larger physical voxel spaces."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry,voxel_pose_terms,SupportedVoxelOperator,VoxelObservationOperator,VoxelReference
from fourier_splats.uq_pose_audit import audit_pose_weights
from fourier_splats.uncertainty import optimize_certificate,bias_aware_half_width
from audit_uq_grid_refinement import model,target
ROOT=Path(__file__).resolve().parents[1];MAPS={'10028':'2660','10049':'6487','10076':'8434'}


def run(args,dataset):
    g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=args.particles,seed=609315)
    ck=np.load(ROOT/'results/uncertainty/development/representation'/dataset/'real_particles-spacing-2.0.npz')
    box=args.coarse_box;op,pilot,mask,noise=model(g,ck,box);B=2.;n=op.op.n;m=2*op.op.q
    ell=target(box,mask,.07,args.target);a=np.deg2rad(args.angle);shift=.01
    terms=voxel_pose_terms(op,pilot,g['q'],a,shift,B,order=2)
    fit=optimize_certificate(op.as_linear_operator(),terms['jacobian'],ell,B,1,terms['remainder'],
        extra_nuisance_groups=[(terms['interaction_factor'],B),(terms['quadratic_jacobian'].reshape(n,m,25),.5),
                               (terms['quadratic_interaction_factor'],.5*B)],gram_diagonal=op.gram_diagonal,maxiter=500,
        rtol=.005,cg_maxiter=500,cg_rtol=1e-7)
    out=ROOT/'results/uncertainty/development'/args.output;out.mkdir(parents=True,exist_ok=True)
    np.savez(out/f'{dataset}-coarse-weights.npz',weights=fit.weights,ell=ell,pilot=pilot,indices=g['indices'])
    result={'dataset':dataset,'config':vars(args),'stage':'development; same fixed weights audited in larger finite density/pose classes',
            'coarse_half_width':fit.half_width,'coarse_width_fraction':float(fit.half_width/(B*np.linalg.norm(ell))),
            'coarse_converged':fit.converged,'coarse_iterations':fit.iterations,'noise_std_physical':float(noise),
            'weights_sha256':hashlib.sha256((out/f'{dataset}-coarse-weights.npz').read_bytes()).hexdigest(),'grids':[]}
    for size in map(int,args.audit_boxes.split(',')):
        begin=time.perf_counter();fine,pilotfine,finemask,_=model(g,ck,size,noise);ellfine=target(size,finemask,.07,args.target)
        audit=audit_pose_weights(fine,pilotfine,g['q'],fit.weights,a,shift,B,order=2,backend=args.backend)
        bias=B*np.linalg.norm(ellfine-fine.adjoint(fit.weights))+sum(audit['totals'].values());sd=np.linalg.norm(fit.weights)
        half=bias_aware_half_width(sd,bias)
        ref=VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map',box=size).volume.ravel()[finemask];ref/=np.linalg.norm(ref)
        u=np.zeros((n,3));u[:,0]=a
        kp=np.einsum('nqi,nij->nqj',g['k'],Rotation.from_rotvec(u).as_matrix())
        pert=SupportedVoxelOperator(VoxelObservationOperator(kp,g['ctf'],size,noise*size**1.5),finemask)
        expected=float(ellfine@pilotfine+fit.weights@(pert.forward(ref)-fine.forward(pilotfine)));truth=float(ellfine@ref)
        error=expected-truth
        cov=lambda q:float(norm.cdf((q-error)/sd)-norm.cdf((-q-error)/sd))
        record={'box':size,'active_voxels':fine.shape[1],'audited_half_width':half,'audited_width_fraction':float(half/(B*np.linalg.norm(ellfine))),
                'width_ratio_to_coarse':float(half/fit.half_width),'bias_components':audit['totals'],
                'density_bias':float(B*np.linalg.norm(ellfine-fine.adjoint(fit.weights))),
                'reference_pilot_distance':float(np.linalg.norm(ref-pilotfine)),'reference_target':truth,
                'reference_analytic_coverage':cov(half),'coarse_width_reference_coverage_with_fine_center':cov(fit.half_width),
                'seconds':time.perf_counter()-begin}
        result['grids'].append(record);(out/f'{dataset}.json').write_text(json.dumps(result,indent=2)+'\n')
        print(dataset,size,'width ratio',record['width_ratio_to_coarse'],'seconds',record['seconds'],flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--particles',type=int,default=128)
    p.add_argument('--coarse-box',type=int,default=24);p.add_argument('--audit-boxes',default='24,32,64');p.add_argument('--angle',type=float,default=.5)
    p.add_argument('--target',choices=['center','contrast'],default='contrast');p.add_argument('--output',default='pose-grid-audit')
    p.add_argument('--backend',choices=['auto','direct','nufft'],default='auto');args=p.parse_args()
    for dataset in args.datasets.split(','):run(args,dataset)
