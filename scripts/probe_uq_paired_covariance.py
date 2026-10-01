#!/usr/bin/env python3
"""Declared finite-orbit matrix cross-covariance feasibility screen."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from fourier_splats.uq_paired_covariance import covariance_cone_diagnostic, direction_growth

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
SOURCES=['src/fourier_splats/uq_paired_covariance.py','tests/test_paired_covariance.py',
    'scripts/probe_uq_paired_covariance.py','research/uncertainty/paired-power-v1/MATRIX-PROPOSAL.md']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise RuntimeError('Commit protocol and source before the molecular screen')
    out=BASE/'paired-covariance-finite-view-v1'
    if out.exists():raise RuntimeError('Preserve every attempt')
    out.mkdir();start=time.perf_counter()
    record=dict(complete=False,cases=[],sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        input_hashes={},arrays={},scope='Oracle finite-view known-translation feasibility screen only.')
    def save():(out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    summary_path=BASE/'paired-power-enlarged-cone-v1/summary.json'
    old=json.loads(summary_path.read_text());record['input_hashes'][str(summary_path.relative_to(ROOT))]=sha(summary_path)
    for ds in ['10028','10049','10076']:
        path=ROOT/old['arrays'][ds]['path'];assert sha(path)==old['arrays'][ds]['sha256']
        record['input_hashes'][str(path.relative_to(ROOT))]=sha(path)
        with np.load(path) as data:
            f=data['fourier'][[0,2]].copy();transfer=np.sqrt(data['transfer_squared'][0])
        means=np.concatenate([f.real,f.imag],axis=-1)*np.tile(transfer,2)
        aggregate=sum(m.T@m/len(m) for m in means)/2
        eigenvalues,vectors=np.linalg.eigh(aggregate);basis=vectors[:,::-1][:,:32]
        compressed=means@basis
        arrays=dict(basis=basis,means=compressed,transfer=transfer,
            basis_eigenvalues=eigenvalues[::-1],full_true_trace=np.array(np.mean(np.sum(means[0]**2,axis=1))))
        for rank in [8,16,32]:
            target=compressed[0,:,:rank];s=target.T@target/len(target)
            for index,name in enumerate(['true_map','region_removed']):
                tick=time.perf_counter();key=f'{name}_{rank}'
                case=dict(dataset=ds,candidate=name,rank=rank,views=len(target),
                    retained_true_energy=float(np.trace(s)/arrays['full_true_trace']))
                try:
                    result=covariance_cone_diagnostic(compressed[index,:,:rank],s,time_limit=60.)
                    for field in ['coefficients','approximation','direction','raw_direction']:
                        arrays[key+'_'+field]=result.pop(field)
                    arrays[key+'_signal_second_moment']=s
                    values={}
                    for variance in [1.,2.,4.]:
                        growth=direction_growth(arrays[key+'_direction'],s,variance)
                        arrays[key+f'_weights_v{variance:g}']=growth.pop('weights')
                        growth['expected_log_128_lower']=128*growth['expected_log_lower']
                        values[str(variance)]=growth
                    case.update(**result,growth=values)
                except Exception as error:
                    case['error']=repr(error)
                case['seconds']=time.perf_counter()-tick;record['cases'].append(case);save()
                print(ds,name,rank,'residual',case.get('scaled_residual'),
                    'growth',case.get('growth',{}).get('1.0',{}).get('expected_log_128_lower'),
                    'error',case.get('error'),'seconds',round(case['seconds'],2),flush=True)
        output=out/f'{ds}-arrays.npz';np.savez_compressed(output,**arrays)
        record['arrays'][ds]=dict(path=str(output.relative_to(ROOT)),sha256=sha(output));save()
    record.update(complete=True,seconds=time.perf_counter()-start);save()
    print('COMPLETE',len(record['cases']),record['seconds'],flush=True)


if __name__=='__main__':main()
