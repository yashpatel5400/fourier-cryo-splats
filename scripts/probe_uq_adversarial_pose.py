#!/usr/bin/env python3
"""Eighty declared local pose searches, each verified by exact cell forward."""
import hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous import cell_forward
from fourier_splats.uq_data import VoxelReference
from fourier_splats.uq_pose_cone_probe import (InterpolatedCellFourier,normalized_quadratic,guide_pose_score)
from fourier_splats.uq_paired_covariance import direction_growth
from fourier_splats.uq_shift_moments import shifted_real_means

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_adversarial_pose.py','src/fourier_splats/uq_pose_cone_probe.py',
 'tests/test_pose_cone_probe.py','research/uncertainty/paired-power-v1/ADVERSARIAL-POSE-PROTOCOL.md']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit source and protocol first')
    out=BASE/'paired-covariance-adversarial-10049-v1'
    if out.exists():raise RuntimeError('Preserve previous attempts')
    out.mkdir();start=time.perf_counter();record=dict(complete=False,cases=[],input_hashes={},
        sources={n:sha(ROOT/n) for n in SOURCES},git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def save():(out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    save();sp=BASE/'paired-covariance-shifted-v1/summary.json';previous=json.loads(sp.read_text())
    ap=ROOT/previous['arrays']['10049']['path'];assert sha(ap)==previous['arrays']['10049']['sha256']
    with np.load(ap) as f:data={k:f[k].copy() for k in f.files}
    rp=BASE/'mixture-validation-preflight-v2/10049.json';geometry=json.loads(rp.read_text())
    map_path=ROOT/'data/uncertainty/references/emd_6487.map'
    for p in [sp,ap,rp,map_path]:record['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
    rho=VoxelReference.from_mrc(map_path,box=64).volume.ravel();rho/=np.linalg.norm(rho)
    x=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij');xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
    center=np.asarray(geometry['target']['center_fraction_field'])
    density=rho*(1-np.exp(-.5*np.sum((xyz-center)**2,axis=1)/(20/geometry['field_A'])**2))
    guide=InterpolatedCellFourier(density,64,256);q=data['q'];transfer=data['transfer'];plane=np.pad(q,((0,0),(0,1)))
    rng=np.random.default_rng(261009);random_rotations=Rotation.random(8,random_state=rng).as_matrix();random_shifts=rng.uniform(-32,32,(8,2))
    archive=dict(random_rotations=random_rotations,random_shifts=random_shifts)
    record['setup_seconds']=time.perf_counter()-start;save()
    for sigma in [0.,.5,1.,2.]:
        key=f'sigma{sigma:g}';direction=data[key+'_union_direction'];target=data[key+'_signal_second_moment']
        means=shifted_real_means(data['fresh_fourier']*transfer,q,sigma*data['fresh_standard_shifts'])
        scores=np.sum((means@direction)*means,axis=1)/np.sum(means*means,axis=1)
        indices=np.argsort(scores)[-12:][::-1]
        starts=[(data['fresh_rotations'][i],sigma*data['fresh_standard_shifts'][i],int(i)) for i in indices]
        starts += [(r,s,None) for r,s in zip(random_rotations,random_shifts)]
        case=dict(shift_sd_pixels=sigma,starts=[],maximum_input_catalog_score=float(scores.max()))
        record['cases'].append(case);save()
        verified_scores=[];selected_rotations=[];selected_shifts=[]
        for j,(base,shift,index) in enumerate(starts):
            tick=time.perf_counter();initial=np.r_[np.zeros(3),shift]
            objective=lambda x:-guide_pose_score(x,base,guide,q,transfer,direction)
            initial_value=-objective(initial)
            fit=minimize(objective,initial,method='L-BFGS-B',bounds=[(-np.pi,np.pi)]*3+[(-32.,32.)]*2,
                options=dict(maxiter=150,maxfun=2000,ftol=1e-12,gtol=1e-8,maxls=30))
            selected=fit.x if -fit.fun>=initial_value else initial
            rotation=base@Rotation.from_rotvec(selected[:3]).as_matrix();k=(plane@rotation)[None]
            values=cell_forward(k,np.ones((1,len(q))),density,64,1.).reshape(1,-1)[0]
            exact=values[:len(q)]+1j*values[len(q):]
            score=normalized_quadratic(exact,q,transfer,direction,selected[3:])
            guided=guide_pose_score(selected,base,guide,q,transfer,direction)
            row=dict(start=j,source_fresh_index=index,status=int(fit.status),success=bool(fit.success),message=str(fit.message),
                evaluations=int(fit.nfev),iterations=int(fit.nit),initial_guide_score=initial_value,
                selected_guide_score=guided,verified_cell_score=score,guide_score_discrepancy=abs(guided-score),
                retained_initial=bool(np.array_equal(selected,initial)),selected_parameters=selected.tolist(),seconds=time.perf_counter()-tick)
            case['starts'].append(row);verified_scores.append(score);selected_rotations.append(rotation);selected_shifts.append(selected[3:]);save()
        correction=max(0.,max(verified_scores))+1e-12;repaired=direction-correction*np.eye(len(direction));after={}
        archive[key+'_selected_rotations']=np.stack(selected_rotations);archive[key+'_selected_shifts']=np.stack(selected_shifts)
        archive[key+'_repaired_direction']=repaired
        for v in [1.,2.,4.]:
            growth=direction_growth(repaired,target,v);archive[key+'_weights_v'+str(v)]=growth.pop('weights');after[str(v)]=growth
        case.update(maximum_verified_violation=max(verified_scores),positive_verified_poses=int(np.sum(np.array(verified_scores)>1e-10)),
            maximum_guide_discrepancy=max(r['guide_score_discrepancy'] for r in case['starts']),
            additional_repair=correction,growth_after_verified_repair=after)
        save();print(sigma,'violation',case['maximum_verified_violation'],'positive',case['positive_verified_poses'],
            'guide_error',case['maximum_guide_discrepancy'],'v1growth',after['1.0']['expected_log_lower'],flush=True)
    path=out/'witnesses.npz';np.savez_compressed(path,**archive)
    record.update(complete=True,seconds=time.perf_counter()-start,arrays=dict(path=str(path.relative_to(ROOT)),sha256=sha(path)));save()
    print('COMPLETE',record['seconds'],flush=True)


if __name__=='__main__':main()
