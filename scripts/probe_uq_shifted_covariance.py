#!/usr/bin/env python3
"""Declared 12-cell common-shift and larger-catalog information gate."""
import hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous import cell_forward
from fourier_splats.uq_data import VoxelReference
from fourier_splats.uq_covariance_cone_bounds import full_covariance_cone,unrestricted_matrix_growth_upper
from fourier_splats.uq_paired_covariance import direction_growth
from fourier_splats.uq_shift_moments import gaussian_shift_second_moment,shifted_real_means

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_shifted_covariance.py','src/fourier_splats/uq_shift_moments.py',
 'src/fourier_splats/uq_covariance_cone_bounds.py','src/fourier_splats/uq_paired_covariance.py',
 'tests/test_shift_moments.py','research/uncertainty/paired-power-v1/SHIFTED-COVARIANCE-PROTOCOL.md']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit protocol and source first')
    out=BASE/'paired-covariance-shifted-v1'
    if out.exists():raise RuntimeError('Preserve every attempt')
    out.mkdir();start=time.perf_counter()
    record=dict(complete=False,cases=[],input_hashes={},arrays={},sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def save():(out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    def remember(path):record['input_hashes'][str(path.relative_to(ROOT))]=sha(path)
    oldp=BASE/'paired-power-enlarged-cone-v1/summary.json';previousp=BASE/'paired-covariance-full-frequency-v1/summary.json'
    old=json.loads(oldp.read_text());previous=json.loads(previousp.read_text());remember(oldp);remember(previousp);save()
    for ds,emd in [('10028','2660'),('10049','6487'),('10076','8434')]:
        arrays={};tick=time.perf_counter()
        try:
            path=ROOT/old['arrays'][ds]['path'];assert sha(path)==old['arrays'][ds]['sha256'];remember(path)
            with np.load(path) as data:
                original=data['fourier'][[0,2]].copy();transfer=np.sqrt(data['transfer_squared'][0]);q=data['q'].copy();oldk=data['original_k'].copy()
            path=ROOT/previous['arrays'][ds]['path'];assert sha(path)==previous['arrays'][ds]['sha256'];remember(path)
            with np.load(path) as data:training=np.concatenate([original[1],data['fresh_fourier']],axis=0)*transfer
            true=original[0]*transfer
            path=BASE/f'mixture-validation-preflight-v2/{ds}.json';remember(path);geometry=json.loads(path.read_text())
            path=ROOT/f'data/uncertainty/references/emd_{emd}.map';remember(path)
            rho=VoxelReference.from_mrc(path,box=64).volume.ravel();rho/=np.linalg.norm(rho)
            x=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij');xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
            center=np.asarray(geometry['target']['center_fraction_field'])
            density=rho*(1-np.exp(-.5*np.sum((xyz-center)**2,axis=1)/(20/geometry['field_A'])**2))
            nq=len(q);f=cell_forward(oldk,np.ones(oldk.shape[:2]),density,64,1.).reshape(len(oldk),2*nq)
            discrepancy=float(np.max(abs(f[:,:nq]+1j*f[:,nq:]-original[1,:64])))
            if discrepancy>1e-8:raise ArithmeticError('Original map construction no longer matches')
            rotations=Rotation.random(10000,random_state=np.random.default_rng(261007+int(ds))).as_matrix()
            fresh=np.empty((len(rotations),nq),complex);plane=np.pad(q,((0,0),(0,1)))
            for begin in range(0,len(rotations),256):
                k=np.einsum('qi,nij->nqj',plane,rotations[begin:begin+256])
                f=cell_forward(k,np.ones(k.shape[:2]),density,64,1.).reshape(len(k),2*nq)
                fresh[begin:begin+len(k)]=f[:,:nq]+1j*f[:,nq:]
            training_shifts=np.random.default_rng(261006+int(ds)).normal(size=(len(training),2))
            fresh_shifts=np.random.default_rng(261008+int(ds)).normal(size=(len(fresh),2))
            arrays.update(fresh_fourier=fresh,fresh_rotations=rotations,training_standard_shifts=training_shifts,
                fresh_standard_shifts=fresh_shifts,q=q,transfer=transfer)
            for sigma in [0.,.5,1.,2.]:
                cell_start=time.perf_counter();key=f'sigma{sigma:g}'
                case=dict(dataset=ds,shift_sd_pixels=sigma,shift_sd_A=sigma*geometry['field_A']/64,
                    training_views=len(training),fresh_views=len(fresh),source_template_discrepancy=discrepancy)
                try:
                    target=gaussian_shift_second_moment(true,q,sigma)
                    null=shifted_real_means(training,q,sigma*training_shifts)
                    fit=full_covariance_cone(null,target,maximum_seconds=120.,maximum_iterations=3000)
                    for v,w in fit.pop('weights').items():arrays[key+'_catalog_weights_v'+v]=w
                    for name in ['coefficients','approximation','direction']:arrays[key+'_'+name]=fit.pop(name)
                    arrays[key+'_signal_second_moment']=target
                    freshmeans=shifted_real_means(fresh*transfer,q,sigma*fresh_shifts)
                    direction=arrays[key+'_direction'];norms=np.sum(freshmeans*freshmeans,axis=1)
                    violations=np.sum((freshmeans@direction)*freshmeans,axis=1)/norms
                    correction=max(0.,float(violations.max()))+1e-12
                    repaired=direction-correction*np.eye(len(direction));after={}
                    for v in [1.,2.,4.]:
                        growth=direction_growth(repaired,target,v)
                        arrays[key+'_union_weights_v'+str(v)]=growth.pop('weights');after[str(v)]=growth
                    arrays[key+'_union_direction']=repaired
                    case.update(fit=fit,maximum_fresh_normalized_violation=float(violations.max()),
                        fraction_fresh_positive_above_1e_minus10=float(np.mean(violations>1e-10)),
                        union_direction_repair=correction,growth_after_union_repair=after,
                        true_map_continuous_mixture_control=unrestricted_matrix_growth_upper(np.zeros_like(target)))
                except Exception as error:case['error']=repr(error)
                case['seconds']=time.perf_counter()-cell_start;record['cases'].append(case);save()
                print(ds,sigma,'catalog',case.get('fit',{}).get('growth',{}).get('1.0'),
                    'fresh_max',case.get('maximum_fresh_normalized_violation'),'after',case.get('growth_after_union_repair',{}).get('1.0'),
                    'error',case.get('error'),'seconds',round(case['seconds'],2),flush=True)
        except Exception as error:
            record.setdefault('setup_errors',[]).append(dict(dataset=ds,error=repr(error),seconds=time.perf_counter()-tick))
            print(ds,'SETUP ERROR',repr(error),flush=True)
        path=out/f'{ds}-arrays.npz';np.savez_compressed(path,**arrays)
        record['arrays'][ds]=dict(path=str(path.relative_to(ROOT)),sha256=sha(path));save()
    record.update(complete=not record.get('setup_errors') and len(record['cases'])==12 and not any('error' in c for c in record['cases']),
        seconds=time.perf_counter()-start);save()
    print('COMPLETE',len(record['cases']),record['seconds'],record['complete'],flush=True)


if __name__=='__main__':main()
