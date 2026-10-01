#!/usr/bin/env python3
"""Frozen Monte Carlo view-law feasibility study; no experimental calibration."""
import hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous import CellObservationOperator
from fourier_splats.uq_data import VoxelReference
from fourier_splats.uq_moment_mc import (noisy_amplitude_polynomial,cubic_interval_maximum,
    binomial_upper,binomial_interval,rejection_count,projected_rejection_probability)
from probe_uq_bispectrum_poses import direct_cells

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_moment_mc.py','scripts/probe_uq_bispectrum_poses.py',
    'src/fourier_splats/uq_moment_mc.py','src/fourier_splats/uq_continuous.py',
    'src/fourier_splats/uq_data.py',
    'research/uncertainty/paired-power-v1/MONTE-CARLO-VIEW-LAW-PROTOCOL.md']
TOTAL=131072;BATCH=512


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit source and protocol first')
    out=BASE/'bispectrum-monte-carlo-v1'
    if out.exists():raise ValueError('Preserve every attempt')
    out.mkdir();start=time.perf_counter()
    result=dict(complete=False,cases=[],stages=[],arrays={},input_hashes={},
        samples_per_stage=TOTAL,batch=BATCH,alpha=.05,calibration_delta=.001,
        sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def save():(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    def read(p):
        result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        return json.loads(p.read_text())
    save();prior=read(BASE/'bispectrum-hull-v1/summary.json')
    for ds,emd in [('10028','2660'),('10049','6487'),('10076','8434')]:
        ap=ROOT/prior['arrays'][ds]['path'];assert sha(ap)==prior['arrays'][ds]['sha256']
        result['input_hashes'][str(ap.relative_to(ROOT))]=sha(ap)
        with np.load(ap) as f:data={k:f[k].copy() for k in f.files}
        definitions=[]
        for case in prior['cases']:
            if case['dataset']!=ds or case['candidate']!='region_removed':continue
            key=case['features']+('_fixed' if len(case['amplitudes'])==1 else '_range09_11')
            direction=data[key+'_direction']
            if not np.linalg.norm(direction):
                result['cases'].append(dict(dataset=ds,key=key,skipped='Constant zero score; rejection probability zero, no projections or random draws'))
                continue
            triads=data['triads'] if case['features']=='power_bispectrum' else np.empty((0,3),int)
            definitions.append(dict(key=key,triads=triads,direction=direction,
                lower=min(case['amplitudes']),upper=max(case['amplitudes']),threshold=case['alternative_mean']))
        if not definitions:save();continue
        mp=ROOT/f'data/uncertainty/references/emd_{emd}.map';result['input_hashes'][str(mp.relative_to(ROOT))]=sha(mp)
        geometry=read(BASE/f'mixture-validation-preflight-v2/{ds}.json')
        rho=VoxelReference.from_mrc(mp,box=64).volume.ravel();rho/=np.linalg.norm(rho)
        x=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij')
        xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
        center=np.array(geometry['target']['center_fraction_field'])
        removed=rho*(1-np.exp(-.5*np.sum((xyz-center)**2,axis=1)/(20/geometry['field_A'])**2))
        plane=np.pad(data['q'],((0,0),(0,1)));nq=len(plane);transfer=data['transfer']
        archive={};states=[]
        for stage in range(2):
            seed=261012+int(ds)+stage*1_000_000;rng=np.random.default_rng(seed)
            stage_record=dict(dataset=ds,stage=stage,seed=seed,initial_rng=rng.bit_generator.state,
                completed=0,operator_checks=[],complete=False)
            result['stages'].append(stage_record);save()
            score_arrays={name:np.empty((len(definitions),TOTAL)) for name in
                ['true_fixed','removed_fixed','true_envelope','removed_envelope']}
            for begin in range(0,TOTAL,BATCH):
                n=min(BATCH,TOTAL-begin);rotations=Rotation.random(n,random_state=rng).as_matrix()
                k=np.einsum('qi,nij->nqj',plane,rotations)
                operator=CellObservationOperator(k,np.ones((n,nq)),64,1.)
                means=[]
                for label,density in [('true',rho),('removed',removed)]:
                    values=operator.forward(density).reshape(n,2*nq)
                    mean=(values[:,:nq]+1j*values[:,nq:])*transfer;means.append(mean)
                    if begin==0:
                        direct=direct_cells(density,k[0])*transfer
                        difference=float(np.max(abs(direct-mean[0])))
                        if difference>1e-8:raise ArithmeticError('Direct physical-cell projection disagrees')
                        stage_record['operator_checks'].append(dict(map=label,maximum_coordinate_difference=difference))
                        archive[f'stage{stage}_{label}_example_mean']=mean[0]
                        archive[f'stage{stage}_{label}_example_direct']=direct
                if begin==0:archive[f'stage{stage}_example_rotation']=rotations[0]
                noise=rng.normal(size=(n,nq))+1j*rng.normal(size=(n,nq))
                for j,d in enumerate(definitions):
                    for label,mean in zip(['true','removed'],means):
                        coefficients=noisy_amplitude_polynomial(mean,noise,d['triads'],d['direction'])
                        fixed=coefficients.sum(axis=1)
                        envelope,_=cubic_interval_maximum(coefficients,d['lower'],d['upper'])
                        if np.any(envelope+1e-11<fixed):raise ArithmeticError('Amplitude interval containing one failed to dominate')
                        score_arrays[label+'_fixed'][j,begin:begin+n]=fixed
                        score_arrays[label+'_envelope'][j,begin:begin+n]=envelope
                stage_record['completed']=begin+n
                if (begin+n)%8192==0:
                    save();print(ds,'stage',stage,'completed',begin+n,'seconds',round(time.perf_counter()-start,2),flush=True)
            stage_record.update(final_rng=rng.bit_generator.state,complete=True);save()
            archive.update({f'stage{stage}_{name}':value for name,value in score_arrays.items()})
        for j,d in enumerate(definitions):
            threshold=d['threshold'];row=dict(dataset=ds,key=d['key'],threshold=threshold,
                amplitude_interval=[d['lower'],d['upper']],counts={},calibration={},projections=[])
            for stage in range(2):
                for name in ['true_fixed','removed_fixed','true_envelope','removed_envelope']:
                    key=f'stage{stage}_{name}';count=int(np.sum(archive[key][j]>threshold))
                    row['counts'][key]=dict(count=count,total=TOTAL,frequency=count/TOTAL,
                        pointwise_95_interval=binomial_interval(count,TOTAL))
            for candidate in ['true','removed']:
                count=row['counts'][f'stage0_{candidate}_envelope']['count']
                upper=binomial_upper(count,TOTAL,.001)
                row['calibration'][candidate]=dict(event_count=count,probability_upper=upper,delta=.001)
                for kappa in [1.,1.01,1.1,2.]:
                    probability=min(1.,kappa*upper)
                    for particles in [1000,10000,100000]:
                        critical=rejection_count(particles,probability,.049)
                        projection=dict(candidate=candidate,kappa=kappa,particles=particles,
                            calibrated_probability_bound=probability,reject_at_count=critical)
                        for label in ['true_fixed','removed_fixed',candidate+'_envelope']:
                            key=f'stage1_{label}';number=row['counts'][key]['count']
                            projection[label]=projected_rejection_probability(number,TOTAL,particles,critical)
                        row['projections'].append(projection)
            result['cases'].append(row);save()
            print(ds,d['key'],'calibrated removed bound',row['calibration']['removed']['probability_upper'],
                'heldout true frequency',row['counts']['stage1_true_fixed']['frequency'],flush=True)
        p=out/f'{ds}-scores.npz';np.savez_compressed(p,**archive)
        result['arrays'][ds]=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size);save()
    result.update(complete=True,seconds=time.perf_counter()-start);save()
    print('COMPLETE',len(result['cases']),result['seconds'],flush=True)


if __name__=='__main__':main()
