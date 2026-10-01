#!/usr/bin/env python3
"""Equal-noise-budget allocation comparison with fresh held-out images."""
import argparse,hashlib,json,subprocess,time
from pathlib import Path
import mrcfile
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous import CellObservationOperator
from fourier_splats.uq_candidate_score import candidate_region,candidate_moment_direction
from fourier_splats.uq_fisher_score import constrained_fisher_direction
from fourier_splats.uq_bispectrum import moment_features
from fourier_splats.uq_view_risk import risk_baselines,paired_variance_decomposition
from fourier_splats.uq_moment_mc import noisy_amplitude_polynomial,rejection_count,projected_rejection_probability,binomial_interval
from fourier_splats.uq_view_variance import grouped_amplitude_envelope,viewing_probability_bounds
from probe_uq_bispectrum_poses import direct_cells
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_replica_allocation.py','src/fourier_splats/uq_fisher_score.py','src/fourier_splats/uq_view_risk.py','src/fourier_splats/uq_candidate_score.py','src/fourier_splats/uq_continuous.py',
 'src/fourier_splats/uq_bispectrum.py','src/fourier_splats/uq_moment_mc.py','src/fourier_splats/uq_view_variance.py',
 'scripts/probe_uq_bispectrum_poses.py','research/uncertainty/paired-power-v1/REPLICA-ALLOCATION-PROTOCOL.md']
AMPS=[.9,1.,1.1];DELETIONS=[0.,.1,.25,.5,1.];KAPPAS=[1.,1.01,1.1,2.,5.];METHODS=['individual_ratio','grouped_ratio','view_variance','cvar_dkw','cvar_split','mean_only','unpaired_variance']


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--dataset',required=True,choices=['10028','10049','10076']);args=parser.parse_args();ds=args.dataset
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():raise ValueError('Commit source and protocol first')
    out=BASE/'candidate-replica-allocation-v1'/ds
    if out.exists():raise ValueError('Preserve every attempt')
    out.mkdir(parents=True);start=time.perf_counter();result=dict(complete=False,dataset=ds,stages=[],cases=[],input_hashes={},
        amplitudes=AMPS,deletions=DELETIONS,calibration_views=8192,heldout_views=131072,replicas=128,cells=16,
        sources={n:sha(ROOT/n) for n in SOURCES},git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def save():(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    def read(p):
        result['input_hashes'][str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    save();prior=read(BASE/'bispectrum-hull-v1/summary.json');ap=ROOT/prior['arrays'][ds]['path'];assert sha(ap)==prior['arrays'][ds]['sha256']
    result['input_hashes'][str(ap.relative_to(ROOT))]=sha(ap)
    with np.load(ap) as f:q=f['q'].copy();triads=f['triads'].copy();transfer=f['transfer'].copy()
    mp=ROOT/f'results/final/{ds}/gaussian-mean.mrc';result['input_hashes'][str(mp.relative_to(ROOT))]=sha(mp)
    with mrcfile.open(mp) as f:
        rho=np.asarray(f.data,float).copy();field=float(f.voxel_size.x)*64
        assert rho.shape==(64,64,64) and abs(float(f.voxel_size.y)*64-field)<1e-4 and abs(float(f.voxel_size.z)*64-field)<1e-4
    norm=float(np.linalg.norm(rho));rho/=norm;mask,region_record=candidate_region(rho,field);region=rho*mask
    result.update(candidate_original_norm=norm,field_A=field,region=region_record);save()
    archive=dict(q=q,triads=triads,transfer=transfer,candidate_density=rho,region_density=region,region_mask=mask)
    definitions=[];nq=len(q);plane=np.pad(q,((0,0),(0,1)))
    old=read(BASE/f'candidate-fisher-score-v1/{ds}/summary.json');assert old['complete']
    fp=ROOT/old['array']['path'];assert sha(fp)==old['array']['sha256'];result['input_hashes'][str(fp.relative_to(ROOT))]=sha(fp)
    with np.load(fp) as f:
        for key,value in [('candidate_density',rho),('region_density',region),('q',q),('triads',triads),('transfer',transfer)]:np.testing.assert_array_equal(f[key],value)
        for c in old['cases']:
            w=f[c['key']+'_direction'].copy();selected=triads if c['feature_family']=='power_bispectrum' else np.empty((0,3),int)
            archive[c['key']+'_direction']=w;definitions.append(dict(key=c['key'],triads=selected,direction=w,threshold=c['design']['threshold']))
            result['cases'].append(dict(key=c['key'],feature_family=c['feature_family'],design=c['design'],baseline_calibration=c['calibration']))
    save()
    for stage,total,batch in [(1,8192,8),(2,131072,512)]:
        seed=261018+int(ds)+stage*1000000;rng=np.random.default_rng(seed)
        record=dict(stage=stage,total=total,batch=batch,seed=seed,initial_rng=rng.bit_generator.state,completed=0,complete=False,operator_checks=[])
        result['stages'].append(record);save()
        if stage==1:
            archive['cal_grouped']=np.empty((4,total,2));archive['cal_individual']=np.empty((4,total,2))
        else:archive['heldout_scores']=np.empty((4,len(DELETIONS),len(AMPS),total))
        for begin in range(0,total,batch):
            n=min(batch,total-begin);rotations=Rotation.random(n,random_state=rng).as_matrix();k=np.einsum('qi,nij->nqj',plane,rotations)
            op=CellObservationOperator(k,np.ones((n,nq)),64,1.);mm=[]
            for label,density in [('candidate',rho),('region',region)]:
                y=op.forward(density.ravel()).reshape(n,2*nq);m=(y[:,:nq]+1j*y[:,nq:])*transfer;mm.append(m)
                if begin==0:
                    direct=direct_cells(density.ravel(),k[0])*transfer;error=float(np.max(abs(m[0]-direct)))
                    if error>1e-8:raise ArithmeticError('Independent cell sum disagrees')
                    record['operator_checks'].append(dict(map=label,maximum_coordinate_difference=error));archive[f'stage{stage}_{label}_first_direct']=direct
            if begin==0:archive[f'stage{stage}_first_rotation']=rotations[0]
            if stage==1:
                shape=(n,2,128,nq);noise=rng.normal(size=shape)+1j*rng.normal(size=shape)
                if begin==0:archive['cal_first_noise']=noise[0]
                m=np.broadcast_to(mm[0][:,None,None,:],shape).reshape(-1,nq);e=noise.reshape(-1,nq)
                for j,d in enumerate(definitions):
                    c=noisy_amplitude_polynomial(m,e,d['triads'],d['direction']).reshape(n,2,128,4)
                    grouped,individual=grouped_amplitude_envelope(c,.9,1.1,16,d['threshold'])
                    archive['cal_grouped'][j,begin:begin+n]=grouped;archive['cal_individual'][j,begin:begin+n]=individual
            else:
                noise=rng.normal(size=(n,nq))+1j*rng.normal(size=(n,nq))
                if begin==0:archive['heldout_first_noise']=noise[0]
                for di,deletion in enumerate(DELETIONS):
                    m=mm[0]-deletion*mm[1]
                    for j,d in enumerate(definitions):
                        c=noisy_amplitude_polynomial(m,noise,d['triads'],d['direction'])
                        for ai,a in enumerate(AMPS):archive['heldout_scores'][j,di,ai,begin:begin+n]=c[:,0]+a*(c[:,1]+a*(c[:,2]+a*c[:,3]))
            record['completed']=begin+n
            if (begin+n)%(2048 if stage==1 else 16384)==0:
                save();print(ds,'stage',stage,'views',begin+n,'seconds',round(time.perf_counter()-start,2),flush=True)
        record.update(complete=True,final_rng=rng.bit_generator.state)
        save()
    for j,d in enumerate(definitions):
        case=result['cases'][j];case['calibration']=viewing_probability_bounds(archive['cal_grouped'][j],archive['cal_individual'][j],.5,KAPPAS)
        case['classical_calibration']=risk_baselines(archive['cal_grouped'][j],.5,KAPPAS)
        case['decomposition']=paired_variance_decomposition(archive['cal_grouped'][j],.5)
        for b,extra in zip(case['calibration']['bounds'],case['classical_calibration']['bounds']):
            assert b['kappa']==extra['kappa'];b.update(extra)
        events=archive['heldout_scores'][j]>d['threshold'];counts=events.sum(axis=-1);case['counts']=counts.tolist();case['projections']=[];case['repeated_groups']=[]
        for replicas,calibration in [(32,case['baseline_calibration']),(128,case['calibration'])]:
            for b in calibration['bounds']:
                for method in METHODS:
                    p=b[method]
                    for n in [1000,10000,100000]:
                        critical=rejection_count(n,p,.049);row=dict(replicas=replicas,calibration_views=32768 if replicas==32 else 8192,kappa=b['kappa'],method=method,particles=n,probability_bound=p,reject_at_count=critical,outcomes=[])
                        for di,deletion in enumerate(DELETIONS):
                            for ai,a in enumerate(AMPS):row['outcomes'].append(dict(deletion=deletion,amplitude=a,**projected_rejection_probability(int(counts[di,ai]),131072,n,critical)))
                        case['projections'].append(row)
                    n=1024;critical=rejection_count(n,p,.049);row=dict(replicas=replicas,calibration_views=32768 if replicas==32 else 8192,kappa=b['kappa'],method=method,particles=n,probability_bound=p,reject_at_count=critical,groups=128,outcomes=[])
                    for di,deletion in enumerate(DELETIONS):
                        for ai,a in enumerate(AMPS):
                            group_counts=events[di,ai].reshape(-1,n).sum(axis=1);rejected=int(np.sum(group_counts>=critical))
                            row['outcomes'].append(dict(deletion=deletion,amplitude=a,event_counts=group_counts.tolist(),rejected=rejected,fraction=rejected/128,pointwise_95_interval=binomial_interval(rejected,128)))
                    case['repeated_groups'].append(row)
        save()
    p=out/'arrays.npz';np.savez_compressed(p,**archive);result.update(complete=True,seconds=time.perf_counter()-start,array=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size));save()
    print('COMPLETE',ds,result['seconds'],flush=True)


if __name__=='__main__':main()
