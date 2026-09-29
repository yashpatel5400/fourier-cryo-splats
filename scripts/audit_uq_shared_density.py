#!/usr/bin/env python3
"""Post-audit fixed nonlinear certificate weights with shared-density bounds."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_pose_audit import pose_derivative_adjoints
from fourier_splats.uq_shared_density import shared_density_bounds
from fourier_splats.uncertainty import bias_aware_half_width
from audit_uq_grid_refinement import model,target
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def run(args,dataset):
    begin=time.perf_counter();folder=BASE/args.source;source=json.loads((folder/f'{dataset}.json').read_text());cfg=source['config']
    saved=np.load(folder/f'{dataset}-coarse-weights.npz');w=saved['weights'];ell=saved['ell'];pilot=saved['pilot'];n=cfg['particles'];box=args.audit_box or cfg['coarse_box']
    g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=n,seed=609315)
    np.testing.assert_array_equal(g['indices'],saved['indices']);ck=np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
    op,recreated,mask,_=model(g,ck,box,source['noise_std_physical'])
    if box==cfg['coarse_box']:np.testing.assert_allclose(pilot,recreated)
    else:pilot=recreated;ell=target(box,mask,.07,cfg['target'])
    first=[];second=[]
    for i in range(n):
        d1,d2=pose_derivative_adjoints(op,g['q'],w,np.deg2rad(cfg['angle']),.01,i,backend='direct')
        first.append(d1);second.append(.5*d2.reshape(len(pilot),25))
    h=op.adjoint(w)-ell;bounds=shared_density_bounds(h,first,second);B=2.
    grid=next(r for r in source['grids'] if r['box']==box);components=grid['bias_components']
    existing_density=grid['density_bias']+components['linear_density']+components['quadratic_density']
    np.testing.assert_allclose(B*bounds['block_triangle'],existing_density,rtol=1e-9)
    other=components['linear_pilot']+components['quadratic_pilot']+components['remainder'];sd=np.linalg.norm(w)
    rows={key:{'density_bias':B*value,'total_bias':B*value+other,'half_width':bias_aware_half_width(sd,B*value+other)} for key,value in bounds.items()}
    for row in rows.values():
        row['width_fraction_of_old']=row['half_width']/grid['audited_half_width']
        row['width_fraction_of_no_data']=float(row['half_width']/(B*np.linalg.norm(ell)))
    stress=BASE/args.stress/f'{dataset}.json';candidates=[]
    if stress.exists() and box==cfg['coarse_box']:
        stressdata=json.loads(stress.read_text());assert stressdata['weights_sha256']==source['weights_sha256']
        for candidate in stressdata['candidates']:
            bias=candidate['absolute_bias'];half=rows['minimum_valid_bound']['half_width']
            candidates.append({'sign':candidate['sign'],'restart':candidate['restart'],'absolute_bias':bias,
                               'bias_fraction_of_new_upper':bias/rows['minimum_valid_bound']['total_bias'],
                               'analytic_coverage':float(norm.cdf((half-bias)/sd)-norm.cdf((-half-bias)/sd))})
            if bias>rows['minimum_valid_bound']['total_bias']+1e-7:raise AssertionError('Feasible bias exceeds shared-density upper bound')
    result={'stage':'development post-audit; same fixed weights; no claim of optimized spectral weights','dataset':dataset,'config':vars(args),
            'source_config':cfg,'weights_sha256':source['weights_sha256'],'bounds':rows,'noise_sd':float(sd),
            'source_json_sha256':hashlib.sha256((folder/f'{dataset}.json').read_bytes()).hexdigest(),
            'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'src/fourier_splats/uq_shared_density.py']},
            'adversarial_candidates':candidates,'seconds':time.perf_counter()-begin}
    out=BASE/args.output;out.mkdir(parents=True,exist_ok=True);(out/f'{dataset}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(dataset,'width ratio',rows['minimum_valid_bound']['width_fraction_of_old'],'seconds',result['seconds'],flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--source',default='pose-grid-center-0.5')
    p.add_argument('--stress',default='nonlinear-stress-center-0.5');p.add_argument('--output',default='shared-density-center-0.5')
    p.add_argument('--audit-box',type=int);args=p.parse_args()
    for dataset in args.datasets.split(','):run(args,dataset)
