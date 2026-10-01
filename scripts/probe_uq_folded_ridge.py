#!/usr/bin/env python3
"""Review-3 six-target folded-width ridge-path development comparison."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_folded_ridge import folded_ridge_search

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['src/fourier_splats/uq_folded_ridge.py','tests/test_folded_ridge.py',
    'scripts/probe_uq_folded_ridge.py','research/uncertainty/FOLDED-RIDGE-REVIEW3-PROTOCOL.md']
TARGETS=[('center',[[0.,0.,0.]],[1.]),('contrast',[[0.,0.,.08],[0.,0.,-.08]],[1.,-1.])]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit protocol and source before molecular outcomes')
    out=BASE/'folded-ridge-review3-v1'
    if out.exists():raise ValueError('Preserve previous attempts')
    out.mkdir();start=time.perf_counter()
    record=dict(complete=False,cases=[],input_hashes={},sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def save():(out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    for ds in ['10028','10049','10076']:
        generator=BASE/f'local-alignment-calibration-v1/{ds}/generators.npz'
        fits_path=BASE/f'refitting-bias-reanalysis-v1/{ds}/true-pose-fits.json'
        fits=json.loads(fits_path.read_text())
        for path in [generator,fits_path]:record['input_hashes'][str(path.relative_to(ROOT))]=sha(path)
        with np.load(generator) as data:
            k,ctf,noise=data['k'],data['ctf'],float(data['noise_std'])
        gram=QuadratureObservationGram(k,ctf,noise,order=40,preconditioner_rank=1024);gram.nthreads=1
        arrays={}
        for target,centers,signs in TARGETS:
            old=fits[target]['fixed_folded'];initial=old['history'][-1]['ridge'];tick=time.perf_counter()
            row=dict(dataset=ds,target=target,initial_ridge=initial,old_sum_weights_folded_width=old['half_width'],
                gaussian_half_widths={str(tau):fits[target][f'gaussian_tau{tau}_fixed']['half_width'] for tau in [2,1]})
            try:
                fitted=folded_ridge_search(gram,centers,signs,.07,2.,initial_ridge=initial,
                    alpha=.05,relative_tolerance=.005,maximum_evaluations=80,maximum_seconds=180.)
                arrays[target+'_weights']=fitted.pop('weights')
                row.update(fit=fitted,width_ratio_to_old=fitted['half_width']/old['half_width'])
            except Exception as error:row['error']=repr(error)
            row['seconds']=time.perf_counter()-tick;record['cases'].append(row);save()
            print(ds,target,'ratio',row.get('width_ratio_to_old'),'evaluations',row.get('fit',{}).get('evaluations'),
                'gap',row.get('fit',{}).get('relative_gap'),'stop',row.get('fit',{}).get('stop_reason'),
                'error',row.get('error'),flush=True)
        path=out/f'{ds}-weights.npz';np.savez_compressed(path,**arrays)
        record.setdefault('arrays',{})[ds]=dict(path=str(path.relative_to(ROOT)),sha256=sha(path));save()
    record.update(complete=True,seconds=time.perf_counter()-start);save()
    print('COMPLETE',len(record['cases']),record['seconds'],flush=True)


if __name__=='__main__':main()
