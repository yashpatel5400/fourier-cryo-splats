#!/usr/bin/env python3
"""Full-frequency two-sided cone gate and fresh-orientation falsification."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous import cell_forward
from fourier_splats.uq_data import VoxelReference
from fourier_splats.uq_paired_covariance import direction_growth
from fourier_splats.uq_covariance_cone_bounds import (full_covariance_cone,
    unrestricted_matrix_growth_upper,log_growth_variance)

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_full_covariance.py','src/fourier_splats/uq_covariance_cone_bounds.py',
 'tests/test_covariance_cone_bounds.py','research/uncertainty/paired-power-v1/FULL-COVARIANCE-PROTOCOL.md']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit protocol and source first')
    out=BASE/'paired-covariance-full-frequency-v1'
    if out.exists():raise RuntimeError('Preserve every attempt')
    out.mkdir();start=time.perf_counter()
    record=dict(complete=False,cases=[],input_hashes={},arrays={},sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def save():(out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    save();oldpath=BASE/'paired-power-enlarged-cone-v1/summary.json';old=json.loads(oldpath.read_text())
    record['input_hashes'][str(oldpath.relative_to(ROOT))]=sha(oldpath)
    for ds,emd in [('10028','2660'),('10049','6487'),('10076','8434')]:
        tick=time.perf_counter();arrays={}
        try:
            ap=ROOT/old['arrays'][ds]['path'];assert sha(ap)==old['arrays'][ds]['sha256']
            rp=BASE/f'mixture-validation-preflight-v2/{ds}.json';original=json.loads(rp.read_text())
            reference_path=ROOT/f'data/uncertainty/references/emd_{emd}.map'
            for path in [ap,rp,reference_path]:record['input_hashes'][str(path.relative_to(ROOT))]=sha(path)
            with np.load(ap) as f:
                ft=f['fourier'][[0,2]].copy();transfer=np.sqrt(f['transfer_squared'][0]);q=f['q'].copy();oldk=f['original_k'].copy()
            keys={tuple(x) for x in q}
            assert len(keys)==len(q) and (0.,0.) not in keys and not any(tuple(-x) in keys for x in q)
            assert np.max(np.abs(q))<32
            means=np.concatenate([ft.real,ft.imag],axis=-1)*np.tile(transfer,2)
            signal=means[0].T@means[0]/len(means[0])
            arrays.update(signal_second_moment=signal,transfer=transfer,q=q)
            direct=np.einsum('ni,nj->ij',means[0],means[0])/len(means[0])
            control=unrestricted_matrix_growth_upper(signal-direct)
            record['cases'].append(dict(dataset=ds,candidate='true_map',nonredundant_frequencies=True,
                direct_mixture_residual=float(np.linalg.norm(signal-direct,'fro')),growth_upper=control,
                known_pose_unit_amplitude_pair_kl=float(np.mean(np.sum((means[0]-means[1])**2,axis=1)))))
            save()
            result=full_covariance_cone(means[1],signal,maximum_seconds=120.,maximum_iterations=3000)
            weights=result.pop('weights')
            for v,w in weights.items():
                arrays['catalog_weights_v'+v]=w
                result['growth'][v]['log_moment_diagnostics']=log_growth_variance(means[0],w,float(v))
            for key in ['coefficients','approximation','direction']:arrays[key]=result.pop(key)
            case=dict(dataset=ds,candidate='region_removed',stage='catalog_complete',**result)
            record['cases'].append(case);save()
            print(ds,'CATALOG','residual',result['relative_residual_frobenius'],'v1',result['growth']['1.0'],flush=True)
            rho=VoxelReference.from_mrc(reference_path,box=64).volume.ravel();rho/=np.linalg.norm(rho)
            x=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij')
            xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
            center=np.asarray(original['target']['center_fraction_field'])
            bump=np.exp(-.5*np.sum((xyz-center)**2,axis=1)/(20/original['field_A'])**2)
            density=rho*(1-bump);nq=len(q)
            values=cell_forward(oldk,np.ones(oldk.shape[:2]),density,64,1.).reshape(len(oldk),2*nq)
            discrepancy=float(np.max(np.abs(values[:,:nq]+1j*values[:,nq:]-ft[1,:64])))
            if discrepancy>1e-8:raise ArithmeticError('Saved Fourier template mismatch')
            rotations=Rotation.random(10000,random_state=np.random.default_rng(261005+int(ds))).as_matrix()
            plane=np.pad(q,((0,0),(0,1)));fresh=np.empty((10000,nq),complex)
            for begin in range(0,len(rotations),256):
                k=np.einsum('qi,nij->nqj',plane,rotations[begin:begin+256])
                f=cell_forward(k,np.ones(k.shape[:2]),density,64,1.).reshape(len(k),2*nq)
                fresh[begin:begin+len(k)]=f[:,:nq]+1j*f[:,nq:]
            freshmeans=np.concatenate([fresh.real,fresh.imag],axis=-1)*np.tile(transfer,2)
            direction=arrays['direction'];norms=np.sum(freshmeans*freshmeans,axis=1)
            violations=np.sum((freshmeans@direction)*freshmeans,axis=1)/norms
            repair=max(0.,float(violations.max()))+1e-12
            corrected=direction-repair*np.eye(len(direction))
            checks=np.sum((freshmeans@corrected)*freshmeans,axis=1)/norms
            if checks.max()>1e-10:raise ArithmeticError('Fresh union repair failed')
            after={}
            for v in [1.,2.,4.]:
                growth=direction_growth(corrected,signal,v);w=growth.pop('weights')
                arrays['union_weights_v'+str(v)]=w
                growth['log_moment_diagnostics']=log_growth_variance(means[0],w,v)
                after[str(v)]=growth
            arrays.update(fresh_fourier=fresh,fresh_rotations=rotations,union_direction=corrected)
            case.update(stage='complete',fresh_views=10000,original_template_discrepancy=discrepancy,
                maximum_fresh_normalized_violation=float(violations.max()),
                fraction_fresh_violations_above_1e_minus10=float(np.mean(violations>1e-10)),
                fresh_union_direction_repair=repair,maximum_repaired_fresh_violation=float(checks.max()),
                growth_after_union_repair=after,total_seconds=time.perf_counter()-tick)
            print(ds,'FRESH','max',case['maximum_fresh_normalized_violation'],'fraction',case['fraction_fresh_violations_above_1e_minus10'],
                'v1lower',after['1.0']['expected_log_lower'],flush=True)
        except Exception as error:
            record.setdefault('errors',[]).append(dict(dataset=ds,error=repr(error),seconds=time.perf_counter()-tick))
            print(ds,'ERROR',repr(error),flush=True)
        path=out/f'{ds}-arrays.npz';np.savez_compressed(path,**arrays)
        record['arrays'][ds]=dict(path=str(path.relative_to(ROOT)),sha256=sha(path));save()
    record.update(complete=not record.get('errors'),seconds=time.perf_counter()-start);save()
    print('COMPLETE',len(record['cases']),record['seconds'],bool(record['complete']),flush=True)


if __name__=='__main__':main()
