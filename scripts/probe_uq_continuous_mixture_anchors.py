#!/usr/bin/env python3
"""Post-outcome continuous-cover anchor diagnostic; no new observations."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from scipy.special import logsumexp
from fourier_splats.uq_mixture_validation import mixture_envelope_upper
from fourier_splats.uq_mixture_anchor import scaled_anchor_upper
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
PROTOCOL=ROOT/'research/uncertainty/CONTINUOUS-MIXTURE-ANCHOR-PROTOCOL.md'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    files=[Path(__file__),PROTOCOL,ROOT/'src/fourier_splats/uq_mixture_anchor.py',
        ROOT/'src/fourier_splats/uq_mixture_validation.py',ROOT/'tests/test_uq_mixture_anchor.py']
    for p in files:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise ValueError('Commit implementation and protocol before evaluation')
    out=BASE/'continuous-mixture-anchors-v1'
    if out.exists():raise RuntimeError('Preserve previous outcomes')
    out.mkdir();start=time.perf_counter();summary=[]
    for ds in ['10028','10049','10076']:
        path=out/f'{ds}.json'
        record=dict(complete=False,dataset=ds,
            source_snapshot=source_snapshot(ROOT,Path(__file__),[str(p.relative_to(ROOT)) for p in files[1:]]),cases=[])
        def save():path.write_text(json.dumps(record,indent=2)+'\n')
        save()
        try:
            source=BASE/f'continuous-mixture-v1/{ds}.json'
            deadline=time.monotonic()+3600
            while True:
                try:
                    ready=source.exists() and json.loads(source.read_text()).get('complete')
                except json.JSONDecodeError:
                    ready=False
                if ready:break
                if time.monotonic()>deadline:raise TimeoutError('One-hour parent waiting budget reached')
                time.sleep(10)
            original=json.loads(source.read_text())
            if not original.get('scientific_run_complete'):raise ValueError('Completed original computation required')
            array_path=ROOT/original['arrays_file']
            if sha(array_path)!=original['arrays_sha256']:raise ValueError('Original cover changed')
            arr=np.load(array_path)
            record['inputs']={str(source.relative_to(ROOT)):sha(source),str(array_path.relative_to(ROOT)):sha(array_path)}
            lower=max(original['fit']['log_likelihood_lower'],original['oracle_feasible_log_likelihood'])
            normalization=original['fit']['normalization'];saved={}
            for prefix,key in [('best','best_cover_log_envelopes'),('final','leaf_log_envelopes')]:
                begin=time.perf_counter();u=arr[key].T
                fitted=mixture_envelope_upper(u,tolerance=.01,max_iterations=500)
                anchors={'row_max':u.max(axis=1),
                    'primal_weights':logsumexp(u+np.log(fitted.weights)[None,:],axis=1),
                    'upper_weights':logsumexp(u+np.log(fitted.upper_anchor_weights)[None,:],axis=1)}
                bounds={k:scaled_anchor_upper(u,v)+normalization for k,v in anchors.items()}
                bounds.update(independent_image=float(u.max(axis=1).sum()+normalization),
                    gaussian_peak=float(normalization),
                    original=float(original['fit']['log_likelihood_upper']))
                if min(bounds.values())<lower-1e-7:raise ArithmeticError('Upper excludes feasible continuous likelihood')
                case=dict(cover=prefix,cells=u.shape[1],seconds=time.perf_counter()-begin,
                    log_likelihood_uppers=bounds,best_upper=min(bounds.values()),
                    continuous_feasible_lower=original['fit']['log_likelihood_lower'],
                    oracle_feasible_lower=original['oracle_feasible_log_likelihood'],
                    gap_to_best_feasible=min(bounds.values())-lower,
                    envelope_relaxation_primal=fitted.envelope_primal+normalization,
                    envelope_unscaled_dual=fitted.log_likelihood_upper+normalization,
                    envelope_gap=fitted.dual_gap,envelope_converged=fitted.converged,
                    envelope_iterations=fitted.iterations)
                record['cases'].append(case);save();print(ds,case,flush=True)
                saved.update({prefix+'_log_anchor_'+k:v for k,v in anchors.items()})
                saved[prefix+'_primal_weights']=fitted.weights
                saved[prefix+'_upper_weights']=fitted.upper_anchor_weights
            output_arrays=out/f'{ds}-arrays.npz';np.savez_compressed(output_arrays,**saved)
            record.update(complete=True,scientific_run_complete=True,
                arrays_file=str(output_arrays.relative_to(ROOT)),arrays_sha256=sha(output_arrays),
                scope='Post-outcome computation on unchanged continuous covers, no new draws or validation ratios.');save()
            summary.append(dict(dataset=ds,cases=record['cases'],result_sha256=sha(path)))
        except Exception as error:
            record.update(complete=True,scientific_run_complete=False,error=repr(error));save();raise
    (out/'summary.json').write_text(json.dumps(dict(complete=True,records=summary,seconds=time.perf_counter()-start),indent=2)+'\n')


if __name__=='__main__':main()
