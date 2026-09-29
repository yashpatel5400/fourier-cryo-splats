#!/usr/bin/env python3
"""Fixed-pose continuum audit with exact cell and continuous Fourier integrals."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry,VoxelReference
from fourier_splats.uq_subspace import matrix_free_certificate
from fourier_splats.uq_continuous import (continuous_residual_norm,cell_target_coefficients,cell_forward,cell_adjoint,CellObservationOperator)
from fourier_splats.uncertainty import bias_aware_half_width
from audit_uq_grid_refinement import model,target
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development';MAPS={'10028':'2660','10049':'6487','10076':'8434'}


def run(args,dataset):
    g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=128,seed=609315);box=24
    ck=np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz');op,pilot,mask,noise=model(g,ck,box)
    pilot_full=op.expand(pilot);pilot_prediction=cell_forward(g['k'],g['ctf'],pilot_full,box,noise)
    fine=args.reference_box;ref=VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map',box=fine).volume.ravel()
    ref/=np.linalg.norm(ref);reference_prediction=cell_forward(g['k'],g['ctf'],ref,fine,noise)
    # Exact overlap of two orthonormal constant-cell bases on the same cube.
    e=np.linspace(-.5,.5,box+1);f=np.linspace(-.5,.5,fine+1)
    overlap=np.sqrt(box*fine)*np.maximum(0,np.minimum(e[1:,None],f[None,1:])-np.maximum(e[:-1,None],f[None,:-1]))
    projection=np.einsum('ia,jb,kc,abc->ijk',overlap,overlap,overlap,ref.reshape((fine,)*3),optimize=True)
    distance=float(np.sqrt(max(0,2-2*pilot_full@projection.ravel())));B=2.
    result={'dataset':dataset,'stage':'development continuous L2 cube, fixed poses/CTFs and known white Gaussian noise',
            'config':vars(args),'density_class':'L2([-1/2,1/2]^3); radius 2 around unit-norm independent constant-cell pilot',
            'generator':'unit-norm deposited-map constant cells; analytic cell Fourier transform, not point samples',
            'noise_std_physical':float(noise),'reference_pilot_L2_distance':distance,'targets':[],
            'code_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'src/fourier_splats/uq_continuous.py']}}
    out=BASE/args.output;out.mkdir(parents=True,exist_ok=True)
    for width in map(float,args.widths.split(',')):
        for name in args.targets.split(','):
            start=time.perf_counter();centers=[[0,0,0]] if name=='center' else [[0,0,.08],[0,0,-.08]];signs=[1] if name=='center' else [1,-1]
            if args.fitting_space=='cube_cells':
                fit_op=CellObservationOperator(g['k'],g['ctf'],args.fitting_box,noise)
                fit_ell=cell_target_coefficients(args.fitting_box,centers,signs,width)
            else:fit_op=op;fit_ell=target(box,mask,width,name)
            fit=matrix_free_certificate(fit_op,fit_ell,B,rtol=.005);w=fit.weights;sd=np.linalg.norm(w)
            audit=continuous_residual_norm(g['k'],g['ctf'],w,noise,centers,signs,width);bias_bound=B*audit['residual_norm']
            half=bias_aware_half_width(sd,bias_bound);truth=float(cell_target_coefficients(fine,centers,signs,width)@ref)
            center=float(cell_target_coefficients(box,centers,signs,width)@pilot_full+w@(reference_prediction-pilot_prediction));bias=center-truth
            analytic=lambda q,b:float(norm.cdf((q-b)/sd)-norm.cdf((-q-b)/sd))
            row={'target':name,'width_fraction_field':width,'fitting_space':args.fitting_space,
                'coarse_optimization_converged':fit.converged,'audit':audit,'noise_sd':float(sd),'half_width':half,
                'width_fraction_of_continuous_no_data':float(half/(B*audit['target_norm'])),
                'reference_true_target':truth,'reference_expected_center':center,'reference_bias':bias,'reference_coverage':analytic(half,bias),
                'continuous_boundary_coverage':analytic(half,bias_bound),'cell_subspaces':[]}
            for cell_box in map(int,args.cell_boxes.split(',')):
                ell=cell_target_coefficients(cell_box,centers,signs,width);projected=cell_adjoint(g['k'],g['ctf'],w,cell_box,noise)-ell
                norm2=float(projected@projected);missing=audit['residual_norm']**2-norm2
                if missing < -1e-7:raise AssertionError('Galerkin residual exceeds continuous norm')
                cell_half=bias_aware_half_width(sd,B*np.sqrt(norm2))
                row['cell_subspaces'].append({'box':cell_box,'projected_residual_norm':float(np.sqrt(norm2)),
                    'unresolved_within_cell_residual_norm':float(np.sqrt(max(0,missing))),
                    'half_width_assuming_cell_subspace':cell_half,'width_fraction_of_continuous':cell_half/half,
                    'continuous_boundary_coverage_using_cell_width':analytic(cell_half,bias_bound),
                    'reference_coverage_using_cell_width':analytic(cell_half,bias)})
            row['seconds']=time.perf_counter()-start;result['targets'].append(row)
            (out/f'{dataset}.json').write_text(json.dumps(result,indent=2)+'\n')
            print(dataset,width,name,'continuous width/no-data',row['width_fraction_of_continuous_no_data'],'seconds',row['seconds'],flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--targets',default='center,contrast')
    p.add_argument('--widths',default='.03,.07');p.add_argument('--cell-boxes',default='16,24,32,48,64');p.add_argument('--reference-box',type=int,default=64)
    p.add_argument('--output',default='continuous-fixed-pose');p.add_argument('--fitting-space',choices=['sphere_points','cube_cells'],default='cube_cells')
    p.add_argument('--fitting-box',type=int,default=24);args=p.parse_args()
    for dataset in args.datasets.split(','):run(args,dataset)
