#!/usr/bin/env python3
"""Oracle shift-invariant third-moment gate, with exact Gaussian variance."""
import hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from fourier_splats.uq_bispectrum import (frequency_triads,moment_features,
    MomentContrast,moment_hull_projection,amplitude_maximum)

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_bispectrum.py','src/fourier_splats/uq_bispectrum.py',
 'tests/test_bispectrum.py','research/uncertainty/paired-power-v1/BISPECTRUM-PROTOCOL.md']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def profile(contrast,means,amplitudes):
    features=moment_features(means,contrast.triads)
    p=features[:,:contrast.n]@contrast.coefficients[:contrast.n]
    b=features[:,contrast.n:]@contrast.coefficients[contrast.n:]
    coefficients=np.zeros((len(means),5));coefficients[:,2]=p;coefficients[:,3]=b
    maxmean,argmean=amplitude_maximum(coefficients,min(amplitudes),max(amplitudes))
    variance=contrast.variance_polynomial(means)
    maxvariance,argvariance=amplitude_maximum(variance,min(amplitudes),max(amplitudes))
    return maxmean,maxvariance,argmean,argvariance


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit protocol and source first')
    out=BASE/'bispectrum-hull-v1'
    if out.exists():raise RuntimeError('Preserve all outcomes')
    out.mkdir();start=time.perf_counter()
    record=dict(complete=False,cases=[],input_hashes={},arrays={},
        sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def save():(out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    op=BASE/'paired-power-enlarged-cone-v1/summary.json';old=json.loads(op.read_text())
    fp=BASE/'paired-covariance-full-frequency-v1/summary.json';freshrecord=json.loads(fp.read_text())
    for p in [op,fp]:record['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
    for ds in ['10028','10049','10076']:
        ap=ROOT/old['arrays'][ds]['path'];bp=ROOT/freshrecord['arrays'][ds]['path']
        assert sha(ap)==old['arrays'][ds]['sha256'] and sha(bp)==freshrecord['arrays'][ds]['sha256']
        for p in [ap,bp]:record['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        with np.load(ap) as a:
            fourier=a['fourier'][[0,2]].copy();transfer=np.sqrt(a['transfer_squared'][0]);q=a['q'].copy()
        with np.load(bp) as a:fresh=a['fresh_fourier'].copy()
        means=fourier*transfer;fresh*=transfer
        triads=frequency_triads(q);archive=dict(q=q,triads=triads,transfer=transfer)
        for name,selected in [('power',np.empty((0,3),int)),('power_bispectrum',triads)]:
            features=moment_features(means,selected);target=features[0].mean(axis=0)
            for ampname,amplitudes in [('fixed',[1.]),('range09_11',[.9,1.,1.1])]:
                key=name+'_'+ampname;tick=time.perf_counter()
                control_weights=np.zeros(len(means[0])*len(amplitudes))
                j=amplitudes.index(1.);control_weights[j*len(means[0]):(j+1)*len(means[0])]=1/len(means[0])
                control_atoms=np.concatenate([moment_features(a*means[0],selected) for a in amplitudes])
                discrepancy=float(np.linalg.norm(control_weights@control_atoms-target))
                if discrepancy>1e-10:raise ArithmeticError('True-map moment witness mismatch')
                record['cases'].append(dict(dataset=ds,candidate='true_map',features=name,amplitudes=amplitudes,
                    direct_mixture_residual=discrepancy,complete=True))
                archive[key+'_true_control_weights']=control_weights
                del control_atoms
                atoms=np.concatenate([moment_features(a*means[1],selected) for a in amplitudes])
                fit=moment_hull_projection(atoms,target)
                for field in ['weights','approximation','direction']:
                    archive[key+'_'+field]=fit.pop(field)
                archive[key+'_target']=target
                direction=archive[key+'_direction'];contrast=MomentContrast(len(q),selected,direction)
                true_expectation,true_variance=contrast.mean_variance(means[0])
                alternative_mean=float(true_expectation.mean())
                alternative_variance=float(true_variance.mean()+true_expectation.var())
                catalog=profile(contrast,means[1],amplitudes);offgrid=profile(contrast,fresh,amplitudes)
                for prefix,values in [('catalog',catalog),('offgrid',offgrid)]:
                    for field,value in zip(['maximum_mean','maximum_variance','mean_amplitude','variance_amplitude'],values):
                        archive[key+'_'+prefix+'_'+field]=value
                catalog_max=float(catalog[0].max());offgrid_max=float(offgrid[0].max())
                combined_max=max(catalog_max,offgrid_max)
                delta=alternative_mean-combined_max
                null_variance=max(float(catalog[1].max()),float(offgrid[1].max()))
                sample_size=((np.sqrt(19*null_variance)+2*np.sqrt(alternative_variance))/delta)**2 if delta>0 else None
                row=dict(dataset=ds,candidate='region_removed',features=name,amplitudes=amplitudes,
                    dimension=len(target),triads=len(selected),**fit,alternative_mean=alternative_mean,
                    alternative_variance=alternative_variance,catalog_amplitude_profiled_maximum=catalog_max,
                    offgrid_amplitude_profiled_maximum=offgrid_max,catalog_profiled_separation=alternative_mean-catalog_max,
                    combined_profiled_separation=delta,maximum_null_noise_variance=null_variance,
                    finite_catalog_cantelli_sufficient_particles=sample_size,total_seconds=time.perf_counter()-tick,complete=True)
                record['cases'].append(row);save()
                print(ds,key,'distance',[row['distance_lower'],row['distance_upper']],row['status'],
                    'profiled separation',delta,'Cantelli n',sample_size,'seconds',row['total_seconds'],flush=True)
                del atoms
        p=out/f'{ds}-arrays.npz';np.savez_compressed(p,**archive)
        record['arrays'][ds]=dict(path=str(p.relative_to(ROOT)),sha256=sha(p));save()
    record.update(complete=True,seconds=time.perf_counter()-start);save()
    print('COMPLETE',len(record['cases']),record['seconds'],flush=True)


if __name__=='__main__':main()
