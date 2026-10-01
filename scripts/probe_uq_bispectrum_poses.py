#!/usr/bin/env python3
"""Continuous pose counterexamples to finite moment-hull directions."""
import hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from fourier_splats.uq_pose_cone_probe import InterpolatedCellFourier
from fourier_splats.uq_bispectrum import moment_features,amplitude_maximum
from fourier_splats.uq_continuous import cell_forward
from fourier_splats.uq_data import VoxelReference

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_bispectrum_poses.py','src/fourier_splats/uq_bispectrum.py',
 'src/fourier_splats/uq_pose_cone_probe.py',
 'research/uncertainty/paired-power-v1/BISPECTRUM-POSE-PROTOCOL.md']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def score(means,triads,direction,lower,upper):
    features=moment_features(np.atleast_2d(means),triads)
    n=np.shape(means)[-1];coefficients=np.zeros((len(features),5))
    coefficients[:,2]=features[:,:n]@direction[:n]
    coefficients[:,3]=features[:,n:]@direction[n:]
    return amplitude_maximum(coefficients,lower,upper)


def direct_cells(density,k,box=64,batch=32):
    x=(np.arange(box)+.5)/box-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij')
    xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
    values=[]
    for begin in range(0,len(k),batch):
        kk=k[begin:begin+batch]
        values.append(np.exp(-2j*np.pi*(kk@xyz.T))@density/box**1.5*np.prod(np.sinc(kk/box),axis=1))
    return np.concatenate(values)


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit source and protocol first')
    out=BASE/'bispectrum-adversarial-pose-v1'
    if out.exists():raise RuntimeError('Preserve every attempt')
    out.mkdir();start=time.perf_counter()
    record=dict(complete=False,cases=[],input_hashes={},arrays={},sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def save():(out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    save();sp=BASE/'bispectrum-hull-v1/summary.json';prior=json.loads(sp.read_text())
    fp=BASE/'paired-covariance-full-frequency-v1/summary.json';freshrecord=json.loads(fp.read_text())
    for p in [sp,fp]:record['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
    for ds,emd in [('10028','2660'),('10049','6487'),('10076','8434')]:
        ap=ROOT/prior['arrays'][ds]['path'];bp=ROOT/freshrecord['arrays'][ds]['path']
        assert sha(ap)==prior['arrays'][ds]['sha256'] and sha(bp)==freshrecord['arrays'][ds]['sha256']
        mp=ROOT/f'data/uncertainty/references/emd_{emd}.map';gp=BASE/f'mixture-validation-preflight-v2/{ds}.json'
        for p in [ap,bp,mp,gp]:record['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        geometry=json.loads(gp.read_text())
        with np.load(ap) as f:data={k:f[k].copy() for k in f.files}
        with np.load(bp) as f:fresh=f['fresh_fourier'].copy();rotations=f['fresh_rotations'].copy()
        q=data['q'];plane=np.pad(q,((0,0),(0,1)));transfer=data['transfer'];nq=len(q)
        rho=VoxelReference.from_mrc(mp,box=64).volume.ravel();rho/=np.linalg.norm(rho)
        x=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij');xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
        center=np.array(geometry['target']['center_fraction_field'])
        density=rho*(1-np.exp(-.5*np.sum((xyz-center)**2,axis=1)/(20/geometry['field_A'])**2))
        guide=None;archive={}
        random=Rotation.random(8,random_state=np.random.default_rng(261011+int(ds))).as_matrix()
        archive['random_rotations']=random
        for oldcase in prior['cases']:
            if oldcase['dataset']!=ds or oldcase['candidate']!='region_removed':continue
            ampname='fixed' if len(oldcase['amplitudes'])==1 else 'range09_11'
            key=oldcase['features']+'_'+ampname;direction=data[key+'_direction']
            if not np.linalg.norm(direction):
                record['cases'].append(dict(dataset=ds,key=key,skipped='Zero direction; no positive separator'));save();continue
            if guide is None:guide=InterpolatedCellFourier(density,64,256)
            triads=data['triads'] if oldcase['features']=='power_bispectrum' else np.empty((0,3),int)
            lower,upper=min(oldcase['amplitudes']),max(oldcase['amplitudes'])
            catalog_scores=score(fresh*transfer,triads,direction,lower,upper)[0]
            indices=np.argsort(catalog_scores)[-12:][::-1]
            starts=[(rotations[i],int(i)) for i in indices]+[(r,None) for r in random]
            case=dict(dataset=ds,key=key,amplitudes=oldcase['amplitudes'],starts=[],alternative_mean=oldcase['alternative_mean'])
            record['cases'].append(case);save();selected_rotations=[];exact_means=[];verified=[]
            for j,(base,index) in enumerate(starts):
                def objective(v):
                    k=plane@(base@Rotation.from_rotvec(v).as_matrix())
                    return -float(score(guide.values(k)*transfer,triads,direction,lower,upper)[0][0])
                initial=np.zeros(3);first=-objective(initial)
                fit=minimize(objective,initial,method='L-BFGS-B',bounds=[(-np.pi,np.pi)]*3,
                    options=dict(maxiter=150,maxfun=2000,ftol=1e-12,gtol=1e-8,maxls=30))
                chosen=fit.x if -fit.fun>=first else initial
                rotation=base@Rotation.from_rotvec(chosen).as_matrix();k=(plane@rotation)[None]
                values=cell_forward(k,np.ones((1,nq)),density,64,1.)
                mean=(values[:nq]+1j*values[nq:])*transfer
                actual,amplitude=score(mean,triads,direction,lower,upper);guided=-objective(chosen)
                case['starts'].append(dict(start=j,source_catalog_index=index,success=bool(fit.success),
                    status=int(fit.status),message=str(fit.message),iterations=int(fit.nit),evaluations=int(fit.nfev),
                    initial_guide_score=first,selected_guide_score=guided,exact_score=float(actual[0]),
                    selected_amplitude=float(amplitude[0]),guide_discrepancy=abs(float(actual[0])-guided)))
                selected_rotations.append(rotation);exact_means.append(mean);verified.append(float(actual[0]));save()
            worst=int(np.argmax(verified));direct=direct_cells(density,plane@selected_rotations[worst])*transfer
            check=float(score(direct,triads,direction,lower,upper)[0][0])
            discrepancy=abs(check-verified[worst])
            if discrepancy>1e-8:raise ArithmeticError('Independent cell sum disagrees')
            maximum=max(verified)
            case.update(maximum_verified_score=maximum,remaining_mean_separation=oldcase['alternative_mean']-maximum,
                scores_above_alternative=int(np.sum(np.array(verified)>oldcase['alternative_mean']+1e-10)),
                maximum_guide_discrepancy=max(r['guide_discrepancy'] for r in case['starts']),
                worst_start=worst,direct_cell_score=check,direct_cell_discrepancy=discrepancy,complete=True)
            archive[key+'_rotations']=np.stack(selected_rotations);archive[key+'_means']=np.stack(exact_means)
            archive[key+'_direct_worst_mean']=direct
            print(ds,key,'separation',case['remaining_mean_separation'],'scores above alternative',case['scores_above_alternative'],
                'direct discrepancy',discrepancy,flush=True);save()
        p=out/f'{ds}-arrays.npz';np.savez_compressed(p,**archive)
        record['arrays'][ds]=dict(path=str(p.relative_to(ROOT)),sha256=sha(p));save()
    record.update(complete=True,seconds=time.perf_counter()-start);save()
    print('COMPLETE',len(record['cases']),record['seconds'],flush=True)


if __name__=='__main__':main()
