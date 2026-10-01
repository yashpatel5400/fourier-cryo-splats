#!/usr/bin/env python3
"""Independent histogramdd/atan2 replay of recorded-pose summary counts."""
import hashlib
import json
import pickle
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development/stack-nuisance-inventory-v1'


def main():
    destination=BASE/'independent-check.json'
    if destination.exists():raise ValueError('Preserve every check')
    raw=BASE/'summary.json';s=json.loads(raw.read_text());assert s['complete']
    checks=[]
    for d in s['datasets']:
        for rel,v in d['inputs'].items():
            assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==v['sha256']
        base=ROOT/f"background/cryodrgn_empiar/empiar{d['dataset']}/inputs"
        a=np.load(next(base.glob('*.cs')),allow_pickle=False)
        r=pickle.load((base/'poses.pkl').open('rb'))[0]
        # Direct ZYZ formulas; unlike the runner, do not call Rotation.as_euler.
        alpha=np.mod(np.arctan2(r[:,1,2],r[:,0,2]),2*np.pi)/(2*np.pi)
        gamma=np.mod(np.arctan2(r[:,2,1],-r[:,2,0]),2*np.pi)/(2*np.pi)
        so3=np.column_stack([alpha,(r[:,2,2]+1)/2,gamma]).astype(float)
        n=r[:,2,:].copy();n[n[:,2]<0]*=-1
        normal=np.column_stack([np.mod(np.arctan2(n[:,1],n[:,0]),2*np.pi)/(2*np.pi),n[:,2]]).astype(float)
        masks={'all_metadata':np.ones(len(a),bool)}
        if 'published_filter' in d['subsets']:
            z=np.zeros(len(a),bool);z[np.asarray(pickle.load((base/'filtered.ind.pkl').open('rb')),int)]=True;masks['published_filter']=z
        for label,mask in masks.items():
            for key,record in d['subsets'][label]['angular'].items():
                family,b=key.split('-');b=int(b)
                x,bins=(so3,[b]*3) if family=='so3' else (normal,[2*b,b])
                for half in [None,0,1]:
                    chosen=mask if half is None else mask&(a['alignments3D/split']==half)
                    c=np.histogramdd(x[chosen],bins=bins,range=[(0,1)]*len(bins))[0].ravel().astype(int)
                    expected=record['counts'] if half is None else record['half_counts'][half]
                    difference=int(np.max(abs(c-expected)))
                    checks.append(dict(dataset=d['dataset'],subset=label,grid=key,half=half,max_count_difference=difference))
                    if difference:raise ArithmeticError(f'Independent binning mismatch: {checks[-1]}')
            fields=d['subsets'][label]['scalar_fields']
            for name,v in fields.items():
                if name.startswith('derived/'):continue
                q=a[name][mask];assert v['finite']==int(np.isfinite(q).sum())
                assert v['minimum']==float(np.min(q[np.isfinite(q)]))
                assert v['maximum']==float(np.max(q[np.isfinite(q)]))
    result=dict(complete=True,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),summary_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),checks=checks,
                limits='Checks count binning and scalar finite/extreme values. Does not certify latent poses, density caps, noise variance, all quantiles or physical amplitude interpretation.')
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print('complete',len(checks),'all bin counts match')


if __name__=='__main__':main()
