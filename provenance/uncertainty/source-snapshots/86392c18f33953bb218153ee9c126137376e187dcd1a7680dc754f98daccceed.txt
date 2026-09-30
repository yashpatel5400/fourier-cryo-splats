#!/usr/bin/env python3
"""Run the declared full-cover refinement on unchanged parent simulations."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from fourier_splats.uq_mixture_refinement import SelectiveCurvatureOrbit, refine_envelope_mixture
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
PROTOCOL=ROOT/'research/uncertainty/CONTINUOUS-MIXTURE-REFINEMENT-PROTOCOL.md'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    files=[Path(__file__),PROTOCOL,
        ROOT/'research/uncertainty/CONTINUOUS-MIXTURE-THEORY.md',
        ROOT/'research/uncertainty/CONTINUOUS-MIXTURE-CURVATURE-THEORY.md',
        *sorted((ROOT/'src/fourier_splats').glob('*.py')),
        ROOT/'tests/test_uq_mixture_refinement.py',ROOT/'tests/test_uq_mixture_fast.py']
    for p in files:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise ValueError('Commit protocol and implementation before outcomes')
    out=BASE/'continuous-mixture-refinement-v1'
    if out.exists():raise RuntimeError('Preserve earlier outcomes')
    out.mkdir();summaries=[];started=time.perf_counter()
    for dataset in ['10028','10049','10076']:
        start=time.perf_counter();path=out/f'{dataset}.json'
        record=dict(complete=False,dataset=dataset,progress=[],
            scope='Post-outcome computation on unchanged synthetic observations; no experimental calibration or structural-test outcome.',
            source_snapshot=source_snapshot(ROOT,Path(__file__),[str(p.relative_to(ROOT)) for p in files[1:]]))
        def save():path.write_text(json.dumps(record,indent=2)+'\n')
        save()
        try:
            input_path=BASE/f'continuous-mixture-v1/{dataset}.json'
            parent_record=json.loads(input_path.read_text())
            assert parent_record['scientific_run_complete']
            arrays_path=ROOT/parent_record['arrays_file']
            assert sha(arrays_path)==parent_record['arrays_sha256']
            data=np.load(arrays_path)
            orbit=SelectiveCurvatureOrbit(data['plane'],data['centers'],data['sigma'],
                data['coefficients'],data['transfer'],data['observed'])
            fields=['leaf_ids','leaf_lower','leaf_upper','leaf_log_envelopes',
                'support_log_kernels','support_angles','best_feasible_weights']
            parent={key:data[key] for key in fields}
            parent['best_feasible_support_count']=parent_record['fit']['best_feasible_support_count']
            record.update(input_hashes={str(p.relative_to(ROOT)):sha(p) for p in [input_path,arrays_path]},
                original_lower=parent_record['fit']['log_likelihood_lower'],
                original_upper=parent_record['fit']['log_likelihood_upper'])
            save()
            def progress(row):
                record['progress'].append(row);save();print(dataset,row,flush=True)
            result=refine_envelope_mixture(orbit,parent,max_splits=32768,batch_size=64,
                mixture_iterations=100,lower_iterations=200,lower_every=4,
                tolerance=1.,wall_seconds=900.,callback=progress)
            destination=out/f'{dataset}-arrays.npz'
            arrays={key:value for key,value in result.items() if isinstance(value,np.ndarray)}
            np.savez_compressed(destination,**arrays)
            scalar={key:value for key,value in result.items() if not isinstance(value,np.ndarray)}
            record.update(complete=True,scientific_run_complete=True,fit=scalar,
                arrays_file=str(destination.relative_to(ROOT)),arrays_sha256=sha(destination),
                seconds=time.perf_counter()-start)
            save();summaries.append(dict(dataset=dataset,scientific_run_complete=True,
                fit=scalar,result_sha256=sha(path)))
        except Exception as error:
            record.update(complete=True,scientific_run_complete=False,error=repr(error),seconds=time.perf_counter()-start)
            save();summaries.append(dict(dataset=dataset,scientific_run_complete=False,
                error=repr(error),result_sha256=sha(path)))
    (out/'summary.json').write_text(json.dumps(dict(complete=True,records=summaries,
        seconds=time.perf_counter()-started,scope='Continuous likelihood computation only.'),indent=2)+'\n')


if __name__=='__main__':main()
