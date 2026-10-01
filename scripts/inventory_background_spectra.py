#!/usr/bin/env python3
"""Descriptive corner spectra with exact adjustment for removed Nyquist modes."""
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from scipy.fft import dctn
from scipy.linalg import solve_triangular

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/uncertainty/development/background-spectrum-inventory-v1'


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def reference():
    h=np.arange(8)[:,None]-np.arange(8)[None,:]
    k=(64*(h==0)-(-1.)**h)/63
    covariance=np.kron(k,k)
    basis=dctn(np.eye(64).reshape(64,8,8),axes=(-2,-1),norm='ortho').reshape(64,64).T[1:]
    transformed=basis@covariance@basis.T
    return basis,covariance,transformed,np.linalg.cholesky(transformed)


def record(x,dct_power,reference_power):
    x=x.reshape(-1,63);mean=x.mean(0);second=x.T@x/len(x);cov=second-np.outer(mean,mean);scale=float(np.trace(second)/63)
    normalized=second/scale;diag=np.diag(normalized);eig=np.linalg.eigvalsh(normalized)
    off=normalized-np.diag(diag);fourth=np.mean(x**4,axis=0)/np.diag(second)**2
    return dict(patches=len(x),reference_whitened_scale=scale,mean_energy_fraction=float(mean@mean/np.trace(second)),
        diagonal_min=float(diag.min()),diagonal_max=float(diag.max()),eigenvalue_min=float(eig.min()),eigenvalue_max=float(eig.max()),
        off_diagonal_frobenius_fraction=float(np.linalg.norm(off)/np.linalg.norm(normalized)),
        marginal_kurtosis_median=float(np.median(fourth)),marginal_kurtosis_max=float(fourth.max()),
        energy_quantiles={str(q):float(np.quantile(np.mean(x*x,axis=1),q)) for q in [.01,.05,.5,.95,.99]},
        dct_power=np.mean(dct_power.reshape(-1,63),axis=0).tolist(),dct_reference_power=reference_power.tolist(),
        whitened_mean=mean.tolist(),second_moment=second.tolist(),centered_covariance=cov.tolist())


def main():
    sources=['scripts/inventory_background_spectra.py','research/uncertainty/BACKGROUND-SPECTRUM-INVENTORY-PROTOCOL.md']
    for name in sources:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():raise ValueError('Commit protocol and source first')
    if OUT.exists():raise ValueError('Preserve inventory attempts')
    OUT.mkdir(parents=True)
    basis,ref,tref,chol=reference()
    # Independent Fourier basis, matching the original 64-pixel output grid.
    freq=np.arange(-31,32);f=np.exp(2j*np.pi*np.arange(8)[:,None]*freq/64)/np.sqrt(63)
    direct=np.kron((f@f.conj().T).real,(f@f.conj().T).real)
    np.testing.assert_allclose(direct,ref,atol=2e-15,rtol=0)
    result=dict(complete=False,sources={p:sha(ROOT/p) for p in sources},git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        reference_maximum_fourier_difference=float(abs(direct-ref).max()),datasets=[])
    def save():(OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    save()
    for ds in ['10028','10049','10076']:
        folder=ROOT/'data'/ds;ip=folder/'images.npy';idxp=folder/'indices.npy';cp=next((ROOT/'background/cryodrgn_empiar'/('empiar'+ds)/'inputs').glob('*.cs'))
        images=np.load(ip,mmap_mode='r');idx=np.load(idxp);cs=np.load(cp,allow_pickle=False)[idx]
        assert images.shape==(8192,64,64)
        patches=np.stack([images[:,y:y+8,x:x+8] for y,x in [(0,0),(0,56),(56,0),(56,56)]],axis=1).astype(float)
        centered=patches-patches.mean(axis=(-2,-1),keepdims=True)
        features=dctn(centered,axes=(-2,-1),norm='ortho').reshape(8192,4,64)[:,:,1:]
        z=solve_triangular(chol,features.reshape(-1,63).T,lower=True).T.reshape(8192,4,63)
        row=dict(dataset=ds,inputs={str(p.relative_to(ROOT)):sha(p) for p in [ip,idxp,cp,folder/'manifest.json']},
            blob_path_count=len(np.unique(cs['blob/path'])),groups=[],checks=[])
        for label,mask in [('all',np.ones(8192,bool)),('source_half0',cs['alignments3D/split']==0),('source_half1',cs['alignments3D/split']==1)]:
            for corner in ['pooled',0,1,2,3]:
                zz=z[mask] if corner=='pooled' else z[mask,int(corner)]
                ff=features[mask] if corner=='pooled' else features[mask,int(corner)]
                row['groups'].append(dict(subset=label,corner=corner,particles=int(mask.sum()),**record(zz,ff*ff,np.diag(tref))))
        # Independent explicit cosine transform (including orthonormal factors).
        x=np.arange(8);u=np.arange(8)[:,None];c=np.sqrt(2/8)*np.cos(np.pi*u*(x+.5)/8);c[0]/=np.sqrt(2)
        for index in [0,4095,8191]:
            for corner in range(4):
                values=(c@centered[index,corner]@c.T).ravel()[1:]
                zz=np.linalg.solve(chol,values)
                error=float(np.max(abs(zz-z[index,corner])))
                if error>1e-12:raise ArithmeticError('Independent patch transform disagrees')
                row['checks'].append(dict(index=index,corner=corner,maximum_error=error))
        path=OUT/(ds+'-moments.npz')
        np.savez_compressed(path,reference_covariance=ref,dct_reference_covariance=tref,source_halves=cs['alignments3D/split'],particle_patch_energy=np.mean(z*z,axis=2),
                            dct_power=np.mean(features*features,axis=1))
        row['arrays']=dict(path=str(path.relative_to(ROOT)),sha256=sha(path),bytes=path.stat().st_size)
        result['datasets'].append(row);save();print(ds,'done',flush=True)
    result['complete']=True;save()


if __name__=='__main__':main()
