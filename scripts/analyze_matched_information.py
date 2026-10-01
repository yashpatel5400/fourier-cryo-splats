#!/usr/bin/env python3
"""Matched known-pose/moment/noise sensitivity ledger using saved held-out draws."""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous import CellObservationOperator
from fourier_splats.uq_bispectrum import MomentContrast,moment_features
from probe_uq_bispectrum_poses import direct_cells
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def mean_envelope(power,cubic):
    p,b=np.asarray(power),np.asarray(cubic)
    values=[.9**2*p+.9**3*b,1.1**2*p+1.1**3*b]
    stationary=np.divide(-2*p,3*b,out=np.ones_like(p),where=b!=0)
    valid=(b!=0)&(stationary>=.9)&(stationary<=1.1)
    v=stationary**2*p+stationary**3*b
    return np.minimum.reduce([values[0],values[1],np.where(valid,v,np.inf)]),np.maximum.reduce([values[0],values[1],np.where(valid,v,-np.inf)])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--dataset',required=True,choices=['10028','10049','10076']);args=parser.parse_args();ds=args.dataset
    sources=['scripts/analyze_matched_information.py','research/uncertainty/MATCHED-INFORMATION-LEDGER-PROTOCOL.md']
    for name in sources:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():raise ValueError('Commit source and protocol first')
    out=BASE/'matched-information-ledger-v1'/ds
    if out.exists():raise ValueError('Preserve every attempt')
    out.mkdir(parents=True);start=time.perf_counter()
    prior=BASE/'candidate-fisher-score-v1'/ds
    old=json.loads((prior/'summary.json').read_text());assert old['complete']
    result=dict(complete=False,dataset=ds,views=8192,input_hashes={'fisher_summary':sha(prior/'summary.json'),'fisher_arrays':sha(prior/'arrays.npz')},
                sources={p:sha(ROOT/p) for p in sources},git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),cases=[])
    def save():(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    save()
    with np.load(prior/'arrays.npz') as f:
        q=f['q'];triads=f['triads'];transfer=f['transfer'];rho=f['candidate_density'];region=f['region_density']
        archived=f['heldout_scores'];directions={c['key']:f[c['key']+'_direction'] for c in old['cases']}
    n=8192;nq=len(q);plane=np.pad(q,((0,0),(0,1)));rng=np.random.default_rng(261017+int(ds)+2000000)
    means=np.empty((2,n,nq),complex);noise=np.empty((n,nq),complex);rotations=np.empty((n,3,3));checks=[]
    for begin in range(0,n,512):
        r=Rotation.random(512,random_state=rng).as_matrix();rotations[begin:begin+512]=r;k=np.einsum('qi,nij->nqj',plane,r)
        op=CellObservationOperator(k,np.ones((512,nq)),64,1.)
        for i,density in enumerate([rho,region]):
            y=op.forward(density.ravel()).reshape(512,2*nq);m=(y[:,:nq]+1j*y[:,nq:])*transfer
            means[i,begin:begin+512]=m
            if begin==0:
                error=float(np.max(abs(m[0]-direct_cells(density.ravel(),k[0])*transfer)));checks.append(error)
                if error>1e-8:raise ArithmeticError('Independent physical mean check failed')
        noise[begin:begin+512]=rng.normal(size=(512,nq))+1j*rng.normal(size=(512,nq))
    m0,mr=means;deletions=[0.,.1,.25,.5,1.]
    norm=np.sum(abs(mr)**2,axis=1)
    result.update(physical_checks=checks,known_pose_region_energy_mean=float(norm.mean()),known_pose_region_energy_standard_error=float(norm.std(ddof=1)/np.sqrt(n)))
    definitions=[]
    for i,c in enumerate(old['cases']):
        key=c['key'];w=directions[key];selected=triads if key.startswith('power_bispectrum') else np.empty((0,3),int)
        definitions.append((key,w,selected))
        for di,deletion in enumerate(deletions):
            scores=moment_features(m0-deletion*mr+noise,selected,noise_variance=1.)@w
            error=float(np.max(abs(scores-archived[i,di,1,:n])))
            if error>1e-8:raise ArithmeticError('Saved score replay failed')
            checks.append(error)
        if len(selected):
            wp=w.copy();wp[nq:]=0;wb=w.copy();wb[:nq]=0
            definitions.extend([(key+'_power_component',wp,selected),(key+'_bispectrum_component',wb,selected)])
    result['maximum_saved_score_difference']=max(checks)
    for key,w,selected in definitions:
        model=MomentContrast(nq,selected,w);nullmean,nullvar=model.mean_variance(m0)
        p=(abs(m0)**2/2)@w[:nq];b=nullmean-p
        global_lo,global_hi=mean_envelope(np.array([p.mean()]),np.array([b.mean()]))
        view_lo,view_hi=mean_envelope(p,b)
        variance=float(nullvar.mean()+nullmean.var(ddof=0));weight_sum=float(w[:nq].sum())
        row=dict(key=key,coefficient_norm=float(np.linalg.norm(w)),power_weight_sum=weight_sum,
                 null_mean=float(nullmean.mean()),noise_variance_mean=float(nullvar.mean()),view_mean_variance=float(nullmean.var()),null_total_variance=variance,
                 global_amplitude_mean_envelope=[float(global_lo[0]),float(global_hi[0])],
                 view_dependent_amplitude_mean_envelope=[float(view_lo.mean()),float(view_hi.mean())],deletions=[])
        for deletion in deletions:
            amean=moment_features(m0-deletion*mr,selected)@w;gap=float(amean.mean()-nullmean.mean())
            row['deletions'].append(dict(deletion=deletion,alternative_mean=float(amean.mean()),signed_mean_gap=gap,
                squared_separation=gap*gap/variance if variance else None,known_pose_twice_kl=float(deletion**2*norm.mean()),
                white_noise_variance_error_matching_gap=abs(gap/weight_sum) if weight_sum else None,
                outside_global_amplitude_mean_envelope=bool(amean.mean()<global_lo[0] or amean.mean()>global_hi[0]),
                outside_view_dependent_mean_envelope=bool(amean.mean()<view_lo.mean() or amean.mean()>view_hi.mean())))
        result['cases'].append(row);save();print(ds,key,'done',flush=True)
    path=out/'replayed-images.npz';np.savez_compressed(path,means=means,noise=noise,rotations=rotations,region_energy=norm)
    result.update(complete=True,seconds=time.perf_counter()-start,replay=dict(path=str(path.relative_to(ROOT)),sha256=sha(path),bytes=path.stat().st_size));save()
    print('COMPLETE',ds,result['seconds'],flush=True)


if __name__=='__main__':main()
