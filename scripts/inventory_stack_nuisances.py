#!/usr/bin/env python3
"""Descriptive recorded-pose/CTF inventory; not latent-view or noise calibration."""
import hashlib
import json
import pickle
import subprocess
import warnings
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/uncertainty/development/stack-nuisance-inventory-v1'


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def summary(x):
    x=np.asarray(x);finite=np.isfinite(x);v=x[finite]
    return dict(count=int(x.size),finite=int(finite.sum()),nonzero=int(np.count_nonzero(v)),
                minimum=float(v.min()) if len(v) else None,maximum=float(v.max()) if len(v) else None,
                mean=float(v.mean()) if len(v) else None,std=float(v.std()) if len(v) else None,
                quantiles={str(q):float(np.quantile(v,q)) for q in [.01,.05,.25,.5,.75,.95,.99]} if len(v) else {})


def counts_record(bins,half,size):
    c=np.bincount(bins,minlength=size);n=len(bins);p=c/n
    halves=[np.bincount(bins[half==h],minlength=size) for h in [0,1]]
    probabilities=[v/v.sum() for v in halves]
    plugin=float(size*np.sum(p*p)-1)
    return dict(particles=n,bins=size,counts=c.tolist(),half_counts=[v.tolist() for v in halves],
                maximum_empirical_ratio=float(size*p.max()),chi_square_plugin=plugin,
                chi_square_iid_unbiased=float(size*np.sum(c*(c-1))/(n*(n-1))-1),
                uniform_plugin_expectation=(size-1)/n,
                cross_half_chi_square=float(size*np.sum(probabilities[0]*probabilities[1])-1),
                half_total_variation=float(.5*np.sum(abs(probabilities[0]-probabilities[1]))),
                half_maximum_ratios=[float(size*v.max()) for v in probabilities],
                empty_fraction=float(np.mean(c==0)))


def association(x,bins,size):
    x=np.asarray(x,float);good=np.isfinite(x);x=x[good];bins=bins[good]
    v=float(np.var(x))
    if v==0:return dict(status='constant_recorded_field',between_bin_variance_fraction=None)
    c=np.bincount(bins,minlength=size);total=np.bincount(bins,weights=x,minlength=size)
    means=np.divide(total,c,out=np.zeros(size),where=c>0)
    between=float(np.sum(c*(means-x.mean())**2)/len(x))
    return dict(status='descriptive_association',between_bin_variance_fraction=between/v,
                finite=int(len(x)),bin_counts=c.tolist(),bin_means=means.tolist())


def main():
    for name in ['scripts/inventory_stack_nuisances.py','research/uncertainty/STACK-NUISANCE-INVENTORY-PROTOCOL.md']:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit source and protocol before execution')
    if OUT.exists():raise ValueError('Preserve every attempted inventory')
    OUT.mkdir(parents=True)
    result=dict(complete=False,git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                source_sha256=sha(Path(__file__)),datasets=[],scope='Recorded metadata only; no true-pose/view/amplitude/noise uncertainty guarantee')
    def save():(OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    save()
    for ds in ['10028','10049','10076']:
        base=ROOT/f'background/cryodrgn_empiar/empiar{ds}/inputs';cs=next(base.glob('*.cs'))
        data=np.load(cs,allow_pickle=False)
        # These are the already archived author metadata, used by the original downloader.
        r=np.asarray(pickle.load((base/'poses.pkl').open('rb'))[0],float)
        if r.shape!=(len(data),3,3):raise ValueError('Pose inventory mismatch')
        convention=np.max(abs(Rotation.from_rotvec(data['alignments3D/pose']).as_matrix().transpose(0,2,1)-r))
        orthogonal=np.max(abs(r@r.transpose(0,2,1)-np.eye(3)))
        if max(convention,orthogonal)>5e-6:raise ValueError('Pose convention/rotation check failed')
        half=data['alignments3D/split'];assert set(np.unique(half))=={0,1}
        masks={'all_metadata':np.ones(len(data),bool)}
        inputs={str(p.relative_to(ROOT)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in [cs,base/'poses.pkl']}
        if (base/'filtered.ind.pkl').exists():
            idx=np.asarray(pickle.load((base/'filtered.ind.pkl').open('rb')),int)
            mask=np.zeros(len(data),bool);mask[idx]=True;masks['published_filter']=mask
            inputs[str((base/'filtered.ind.pkl').relative_to(ROOT))]=dict(sha256=sha(base/'filtered.ind.pkl'),bytes=(base/'filtered.ind.pkl').stat().st_size)
        with warnings.catch_warnings(record=True) as caught:
            angles=Rotation.from_matrix(r).as_euler('ZYZ')
        angular_warnings=[str(w.message) for w in caught]
        unit=np.column_stack([(angles[:,0]%(2*np.pi))/(2*np.pi),(np.cos(angles[:,1])+1)/2,(angles[:,2]%(2*np.pi))/(2*np.pi)])
        normal=r[:,2,:].copy();normal[normal[:,2]<0]*=-1
        phi=(np.arctan2(normal[:,1],normal[:,0])%(2*np.pi))/(2*np.pi)
        grids={}
        for b in [4,8,12]:
            index=np.clip((unit*b).astype(int),0,b-1)
            grids[f'so3-{b}']=(index[:,0]*b*b+index[:,1]*b+index[:,2],b**3)
            ix=np.clip((phi*2*b).astype(int),0,2*b-1);iz=np.clip((normal[:,2]*b).astype(int),0,b-1)
            grids[f'normal-{b}']=(ix*b+iz,2*b*b)
        row=dict(dataset=ds,inputs=inputs,particles=len(data),rotation_transpose_error=float(convention),
                 orthogonality_error=float(orthogonal),euler_warnings=angular_warnings,subsets={})
        for label,mask in masks.items():
            a=data[mask]
            fields={name:summary(a[name]) for name in data.dtype.names
                    if name.startswith(('ctf/','alignments3D/')) and a[name].ndim==1 and a[name].dtype.kind in 'ifub'}
            fields['derived/defocus_mean_A']=summary((a['ctf/df1_A']+a['ctf/df2_A'])/2)
            fields['derived/astigmatism_abs_A']=summary(abs(a['ctf/df1_A']-a['ctf/df2_A']))
            row['subsets'][label]=dict(particles=int(mask.sum()),scalar_fields=fields,
                angular={name:counts_record(index[mask],half[mask],size) for name,(index,size) in grids.items()},
                association={name:association(a[name],grids['normal-8'][0][mask],128)
                             for name in ['ctf/scale','alignments3D/alpha','alignments3D/weight']})
        result['datasets'].append(row);save()
        print(ds,'complete',len(data),'convention error',convention,flush=True)
    result['complete']=True;save()


if __name__=='__main__':main()
