#!/usr/bin/env python3
"""Fresh known-density preferred-view null and alternative controls."""
import argparse,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from fourier_splats.uq_continuous import CellObservationOperator
from fourier_splats.uq_data import VoxelReference
from fourier_splats.uq_bispectrum import moment_features
from fourier_splats.uq_moment_mc import noisy_amplitude_polynomial,binomial_interval,rejection_count,projected_rejection_probability
from fourier_splats.uq_preferred_views import sample_cap_views
from probe_uq_bispectrum_poses import direct_cells

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_preferred_views.py','src/fourier_splats/uq_preferred_views.py',
 'src/fourier_splats/uq_moment_mc.py','src/fourier_splats/uq_bispectrum.py',
 'src/fourier_splats/uq_continuous.py','src/fourier_splats/uq_data.py',
 'scripts/probe_uq_bispectrum_poses.py','research/uncertainty/paired-power-v1/PREFERRED-VIEW-PROTOCOL.md']
TOTAL=65536;BATCH=512;AMPLITUDES=[.9,1.,1.1];METHODS=['individual_ratio','grouped_ratio','view_variance']


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--dataset',required=True,choices=['10028','10049']);args=parser.parse_args();ds=args.dataset
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():raise ValueError('Commit source and protocol first')
    out=BASE/'bispectrum-preferred-view-v1'/ds
    if out.exists():raise ValueError('Preserve every attempt')
    out.mkdir(parents=True);start=time.perf_counter()
    result=dict(complete=False,dataset=ds,laws=[],input_hashes={},amplitudes=AMPLITUDES,total_per_law=TOTAL,batch=BATCH,
        sources={n:sha(ROOT/n) for n in SOURCES},git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def save():(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    def read(p):
        result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        return json.loads(p.read_text())
    save();prior=read(BASE/'bispectrum-hull-v1/summary.json')
    calibration=read(BASE/f'bispectrum-view-variance-v1/{ds}/summary.json');assert calibration['complete']
    ap=ROOT/prior['arrays'][ds]['path'];assert sha(ap)==prior['arrays'][ds]['sha256'];result['input_hashes'][str(ap.relative_to(ROOT))]=sha(ap)
    with np.load(ap) as f:data={k:f[k].copy() for k in f.files}
    definitions=[]
    for c in calibration['cases']:
        definitions.append(dict(key=c['key'],triads=data['triads'] if c['key'].startswith('power_bispectrum') else np.empty((0,3),int),
            direction=data[c['key']+'_direction'],threshold=c['threshold'],calibration=c['calibration']))
    mp=ROOT/f"data/uncertainty/references/emd_{dict(zip(['10028','10049'],['2660','6487']))[ds]}.map"
    result['input_hashes'][str(mp.relative_to(ROOT))]=sha(mp)
    geometry=read(BASE/f'mixture-validation-preflight-v2/{ds}.json')
    rho=VoxelReference.from_mrc(mp,box=64).volume.ravel();rho/=np.linalg.norm(rho)
    x=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij');xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
    center=np.array(geometry['target']['center_fraction_field'])
    removed=rho*(1-np.exp(-.5*np.sum((xyz-center)**2,axis=1)/(20/geometry['field_A'])**2))
    plane=np.pad(data['q'],((0,0),(0,1)));nq=len(plane);transfer=data['transfer']
    laws=[(axis,kappa) for kappa in [1.1,2.,5.] for axis in range(3)]
    for law_index,(axis,kappa) in enumerate(laws):
        seed=261014+int(ds)+100000*law_index;rng=np.random.default_rng(seed);cut=1-1/kappa
        record=dict(index=law_index,axis=axis,kappa=kappa,seed=seed,cut=cut,complete=False,completed=0,
            proposed=0,accepted=0,used=0,minimum_absolute_coordinate=1.,coordinate_square_sum=0.,
            initial_rng=rng.bit_generator.state,operator_checks=[],moment_checks=[],cases=[])
        result['laws'].append(record);save()
        archive={label:np.empty((len(definitions),len(AMPLITUDES),TOTAL)) for label in ['true','removed']}
        for begin in range(0,TOTAL,BATCH):
            n=min(BATCH,TOTAL-begin);rotations,counts=sample_cap_views(n,axis,kappa,rng)
            for name in ['proposed','accepted','used']:record[name]+=counts[name]
            coordinate=rotations[:,2,axis]
            record['minimum_absolute_coordinate']=min(record['minimum_absolute_coordinate'],float(np.min(abs(coordinate))))
            record['coordinate_square_sum']+=float(coordinate@coordinate)
            k=np.einsum('qi,nij->nqj',plane,rotations);operator=CellObservationOperator(k,np.ones((n,nq)),64,1.)
            noise=rng.normal(size=(n,nq))+1j*rng.normal(size=(n,nq))
            if begin==0:archive.update(example_rotation=rotations[0],example_noise=noise[0])
            for label,density in [('true',rho),('removed',removed)]:
                values=operator.forward(density).reshape(n,2*nq);mean=(values[:,:nq]+1j*values[:,nq:])*transfer
                if begin==0:
                    direct=direct_cells(density,k[0])*transfer;difference=float(np.max(abs(direct-mean[0])))
                    if difference>1e-8:raise ArithmeticError('Physical-cell sum mismatch')
                    record['operator_checks'].append(dict(map=label,maximum_coordinate_difference=difference));archive[label+'_example_direct']=direct
                for j,d in enumerate(definitions):
                    c=noisy_amplitude_polynomial(mean,noise,d['triads'],d['direction'])
                    for a_index,a in enumerate(AMPLITUDES):
                        scores=c[:,0]+a*(c[:,1]+a*(c[:,2]+a*c[:,3]));archive[label][j,a_index,begin:begin+n]=scores
                        if begin==0:
                            direct_score=float(moment_features(a*direct+noise[0],d['triads'],1.)@d['direction'])
                            discrepancy=abs(scores[0]-direct_score)
                            if discrepancy>1e-8:raise ArithmeticError('First direct score mismatch')
                            record['moment_checks'].append(dict(map=label,key=d['key'],amplitude=a,absolute_error=float(discrepancy)))
            record['completed']=begin+n
            if (begin+n)%16384==0:
                save();print(ds,'law',law_index,'views',begin+n,'seconds',round(time.perf_counter()-start,2),flush=True)
        for j,d in enumerate(definitions):
            for a_index,a in enumerate(AMPLITUDES):
                events={label:archive[label][j,a_index]>d['threshold'] for label in ['true','removed']}
                counts={label:int(e.sum()) for label,e in events.items()}
                row=dict(key=d['key'],amplitude=a,threshold=d['threshold'],counts=counts,projections=[],repeated_groups=[])
                for candidate in ['true','removed']:
                    bounds=next(b for b in d['calibration'][candidate]['bounds'] if b['kappa']==kappa)
                    for method in METHODS:
                        p=bounds[method]
                        for n in [1000,10000,100000]:
                            critical=rejection_count(n,p,.049)
                            row['projections'].append(dict(candidate=candidate,method=method,particles=n,probability_bound=p,
                                reject_at_count=critical,**{label:projected_rejection_probability(counts[label],TOTAL,n,critical) for label in ['true','removed']}))
                        n=1024;critical=rejection_count(n,p,.049);replicate=dict(candidate=candidate,method=method,particles=n,groups=TOTAL//n,
                            probability_bound=p,reject_at_count=critical)
                        for label,e in events.items():
                            group_counts=e.reshape(-1,n).sum(axis=1);rejected=int(np.sum(group_counts>=critical))
                            replicate[label]=dict(event_counts=group_counts.tolist(),rejected=rejected,fraction=rejected/len(group_counts),
                                pointwise_95_interval=binomial_interval(rejected,len(group_counts)))
                        row['repeated_groups'].append(replicate)
                record['cases'].append(row)
        p=out/f'law{law_index}-scores.npz';np.savez_compressed(p,**archive)
        record.update(complete=True,final_rng=rng.bit_generator.state,coordinate_square_mean=record['coordinate_square_sum']/TOTAL,
            theoretical_coordinate_square_mean=(1+cut+cut*cut)/3,proposal_acceptance_fraction=record['accepted']/record['proposed'],
            array=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size));save()
        print(ds,'law',law_index,'complete','seconds',round(time.perf_counter()-start,2),flush=True)
    result.update(complete=True,seconds=time.perf_counter()-start);save();print('COMPLETE',ds,result['seconds'],flush=True)


if __name__=='__main__':main()
