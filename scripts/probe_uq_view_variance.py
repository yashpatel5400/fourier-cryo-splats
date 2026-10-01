#!/usr/bin/env python3
"""Noise-replicated view-law follow-up; one independent stack per process."""
import argparse,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous import CellObservationOperator
from fourier_splats.uq_data import VoxelReference
from fourier_splats.uq_bispectrum import moment_features
from fourier_splats.uq_moment_mc import (noisy_amplitude_polynomial,binomial_interval,
    rejection_count,projected_rejection_probability)
from fourier_splats.uq_view_variance import grouped_amplitude_envelope,viewing_probability_bounds
from probe_uq_bispectrum_poses import direct_cells

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_view_variance.py','scripts/probe_uq_bispectrum_poses.py',
    'src/fourier_splats/uq_view_variance.py','src/fourier_splats/uq_moment_mc.py',
    'src/fourier_splats/uq_bispectrum.py','src/fourier_splats/uq_continuous.py',
    'src/fourier_splats/uq_data.py','research/uncertainty/paired-power-v1/VIEW-VARIANCE-PROTOCOL.md']
CALIBRATION=32768;HELDOUT=131072;REPLICAS=32;CELLS=16;KAPPAS=[1.,1.01,1.1,2.,5.]


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--dataset',required=True,choices=['10028','10049']);args=parser.parse_args()
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit source and protocol first')
    ds=args.dataset;out=BASE/'bispectrum-view-variance-v1'/ds
    if out.exists():raise ValueError('Preserve every attempt')
    out.mkdir(parents=True);start=time.perf_counter()
    result=dict(complete=False,dataset=ds,cases=[],stages=[],input_hashes={},
        calibration_views=CALIBRATION,noise_groups=2,replicas_per_group=REPLICAS,
        amplitude_cells=CELLS,heldout_views=HELDOUT,alpha=.05,calibration_delta=.001,
        sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def save():(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    def read(p):
        result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        return json.loads(p.read_text())
    save();prior=read(BASE/'bispectrum-hull-v1/summary.json')
    old_mc=read(BASE/'bispectrum-monte-carlo-v1/summary.json')
    ap=ROOT/prior['arrays'][ds]['path'];assert sha(ap)==prior['arrays'][ds]['sha256']
    result['input_hashes'][str(ap.relative_to(ROOT))]=sha(ap)
    with np.load(ap) as f:data={k:f[k].copy() for k in f.files}
    definitions=[]
    for case in old_mc['cases']:
        if case['dataset']!=ds or not case['key'].endswith('_range09_11'):continue
        key=case['key'];direction=data[key+'_direction'];assert np.linalg.norm(direction)>0
        triads=data['triads'] if key.startswith('power_bispectrum') else np.empty((0,3),int)
        definitions.append(dict(key=key,triads=triads,direction=direction,threshold=case['threshold'],
            centers={label:case['counts'][f'stage0_{label}_envelope']['frequency'] for label in ['true','removed']}))
    assert len(definitions)==2
    emd={'10028':'2660','10049':'6487'}[ds]
    mp=ROOT/f'data/uncertainty/references/emd_{emd}.map';result['input_hashes'][str(mp.relative_to(ROOT))]=sha(mp)
    geometry=read(BASE/f'mixture-validation-preflight-v2/{ds}.json')
    rho=VoxelReference.from_mrc(mp,box=64).volume.ravel();rho/=np.linalg.norm(rho)
    x=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij')
    xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
    center=np.array(geometry['target']['center_fraction_field'])
    removed=rho*(1-np.exp(-.5*np.sum((xyz-center)**2,axis=1)/(20/geometry['field_A'])**2))
    plane=np.pad(data['q'],((0,0),(0,1)));nq=len(plane);transfer=data['transfer'];archive={}
    for stage in range(2):
        total=CALIBRATION if stage==0 else HELDOUT;batch=32 if stage==0 else 512
        seed=261013+int(ds)+stage*1_000_000;rng=np.random.default_rng(seed)
        record=dict(stage=stage,seed=seed,batch=batch,total_views=total,initial_rng=rng.bit_generator.state,
            completed=0,operator_checks=[],complete=False);result['stages'].append(record);save()
        if stage==0:
            for label in ['true','removed']:
                for kind in ['grouped','individual']:archive[f'cal_{label}_{kind}']=np.empty((2,total,2))
        else:
            for label in ['true','removed']:archive[f'heldout_{label}']=np.empty((2,total))
        for begin in range(0,total,batch):
            n=min(batch,total-begin);rotations=Rotation.random(n,random_state=rng).as_matrix()
            k=np.einsum('qi,nij->nqj',plane,rotations)
            operator=CellObservationOperator(k,np.ones((n,nq)),64,1.)
            noise_shape=(n,2,REPLICAS,nq) if stage==0 else (n,nq)
            noise=rng.normal(size=noise_shape)+1j*rng.normal(size=noise_shape)
            for label,density in [('true',rho),('removed',removed)]:
                values=operator.forward(density).reshape(n,2*nq)
                mean=(values[:,:nq]+1j*values[:,nq:])*transfer
                if begin==0:
                    direct=direct_cells(density,k[0])*transfer;difference=float(np.max(abs(direct-mean[0])))
                    if difference>1e-8:raise ArithmeticError('Independent physical-cell sum disagrees')
                    record['operator_checks'].append(dict(map=label,maximum_coordinate_difference=difference))
                    archive[f'stage{stage}_{label}_example_direct']=direct
                    archive[f'stage{stage}_{label}_example_mean']=mean[0]
                if stage==0:
                    repeated=np.broadcast_to(mean[:,None,None,:],noise_shape).reshape(-1,nq)
                    flattened=noise.reshape(-1,nq)
                    for j,d in enumerate(definitions):
                        c=noisy_amplitude_polynomial(repeated,flattened,d['triads'],d['direction']).reshape(n,2,REPLICAS,4)
                        grouped,individual=grouped_amplitude_envelope(c,.9,1.1,CELLS,d['threshold'])
                        archive[f'cal_{label}_grouped'][j,begin:begin+n]=grouped
                        archive[f'cal_{label}_individual'][j,begin:begin+n]=individual
                else:
                    for j,d in enumerate(definitions):
                        archive[f'heldout_{label}'][j,begin:begin+n]=moment_features(mean+noise,d['triads'],1.)@d['direction']
            if begin==0:archive[f'stage{stage}_example_rotation']=rotations[0]
            record['completed']=begin+n
            if (begin+n)%(2048 if stage==0 else 16384)==0:
                save();print(ds,'stage',stage,'views',begin+n,'seconds',round(time.perf_counter()-start,2),flush=True)
        record.update(complete=True,final_rng=rng.bit_generator.state);save()
    for j,d in enumerate(definitions):
        row=dict(key=d['key'],threshold=d['threshold'],amplitude_interval=[.9,1.1],
            calibration={},heldout={},projections=[])
        for label in ['true','removed']:
            count=int(np.sum(archive[f'heldout_{label}'][j]>d['threshold']))
            row['heldout'][label]=dict(count=count,total=HELDOUT,frequency=count/HELDOUT,
                pointwise_95_interval=binomial_interval(count,HELDOUT))
        for candidate in ['true','removed']:
            bound=viewing_probability_bounds(archive[f'cal_{candidate}_grouped'][j],
                archive[f'cal_{candidate}_individual'][j],d['centers'][candidate],KAPPAS)
            row['calibration'][candidate]=bound
            for b in bound['bounds']:
                for method in ['individual_ratio','grouped_ratio','view_variance']:
                    for particles in [1000,10000,100000]:
                        critical=rejection_count(particles,b[method],.049)
                        projection=dict(candidate=candidate,kappa=b['kappa'],method=method,particles=particles,
                            probability_bound=b[method],reject_at_count=critical)
                        for label in ['true','removed']:
                            projection[label]=projected_rejection_probability(row['heldout'][label]['count'],HELDOUT,particles,critical)
                        row['projections'].append(projection)
        result['cases'].append(row);save()
        print(ds,d['key'],'removed calibration',json.dumps(row['calibration']['removed']),flush=True)
    p=out/'arrays.npz';np.savez_compressed(p,**archive)
    result.update(complete=True,seconds=time.perf_counter()-start,
        array=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size));save()
    print('COMPLETE',ds,result['seconds'],flush=True)


if __name__=='__main__':main()
