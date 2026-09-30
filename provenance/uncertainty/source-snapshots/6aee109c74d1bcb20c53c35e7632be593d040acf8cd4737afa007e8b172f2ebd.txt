#!/usr/bin/env python3
"""Post-outcome shared-scale follow-up on the complete discrete-view screen."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from fourier_splats.uq_mixture_scale import profile_common_scale
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
PROTOCOL = ROOT/'research/uncertainty/MIXTURE-COMMON-SCALE-PROTOCOL.md'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    files = [Path(__file__), PROTOCOL,
        ROOT/'src/fourier_splats/uq_mixture_scale.py',
        ROOT/'src/fourier_splats/uq_mixture_validation.py',
        ROOT/'tests/test_uq_mixture_scale.py',
        ROOT/'research/uncertainty/MIXTURE-VALIDATION-PREFLIGHT-RESULTS.md']
    for p in files:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT) != p.read_bytes():
            raise ValueError('Commit protocol and implementation before outcomes')
    out = BASE/'mixture-common-scale-v1'
    if out.exists():
        raise RuntimeError('Preserve previous outcomes')
    out.mkdir(); started = time.perf_counter(); summaries=[]
    for ds in ['10028','10049','10076']:
        path = out/f'{ds}.json'
        record = dict(complete=False,dataset=ds,
            scope='Post-outcome common-scale development on the unchanged discrete-view simulations; not SO(3) or experimental calibration.',
            source_snapshot=source_snapshot(ROOT,Path(__file__),[str(p.relative_to(ROOT)) for p in files[1:]]),cases=[])
        def save():
            path.write_text(json.dumps(record,indent=2)+'\n')
        save()
        try:
            source = BASE/f'mixture-validation-preflight-v2/{ds}.json'
            original = json.loads(source.read_text())
            if not original['complete'] or not original['scientific_run_complete']:
                raise ValueError('Complete original screen required')
            array_path = ROOT/original['arrays_file']
            if sha(array_path) != original['arrays_sha256']:
                raise ValueError('Original simulation arrays changed')
            arrays = np.load(array_path); means = arrays['means']; observations = arrays['observations']
            names = ['true_map','region_half_removed','region_removed','zero_signal']
            record.update(input_hashes={str(source.relative_to(ROOT)):sha(source),str(array_path.relative_to(ROOT)):sha(array_path)},
                dimension=original['dimension'],prior_repeats=len(original['repeats']))
            for rep in original['repeats']:
                y = observations[rep['repeat']]
                for m,name in enumerate(names):
                    if time.perf_counter()-started > 1800:
                        raise TimeoutError('Declared 30-minute budget reached between cases')
                    begin = time.perf_counter()
                    residual = np.sqrt(np.sum((y[:,None,:]-means[m])**2,axis=-1))
                    fitted = profile_common_scale(residual,original['dimension'],
                        tolerance=1.,max_nodes=65,mixture_iterations=2000,mixture_tolerance=1e-4)
                    score = rep['log_oracle_numerator']-fitted['log_likelihood_upper']
                    known = next(c for c in rep['cases'] if c['candidate']==name and c['noise_scope']=='known')
                    separate = next(c for c in rep['cases'] if c['candidate']==name and c['noise_scope']=='profiled_per_image_cell')
                    known_score = known['log_evalues']['shared_law_upper']
                    if score > known_score+known['gap']+2e-6:
                        raise ArithmeticError('Shared-scale upper below a feasible known-scale likelihood')
                    if name == 'true_map' and score > 2e-6:
                        raise ArithmeticError('Matched oracle-numerator control exceeds one')
                    record['cases'].append(dict(repeat=rep['repeat'],candidate=name,
                        log_evalue=float(score),rejected_at_005=bool(score>=np.log(20)),
                        known_scale_log_evalue=known_score,
                        separate_scale_log_evalue=separate['log_evalues']['shared_law_upper'],
                        seconds=time.perf_counter()-begin,fit=fitted))
                    save();print(ds,rep['repeat'],name,round(score,4),round(fitted['gap'],5),len(fitted['nodes']),fitted['status'],flush=True)
            summary=[]
            for name in names:
                cases=[c for c in record['cases'] if c['candidate']==name]
                summary.append(dict(candidate=name,repeats=len(cases),
                    rejections=sum(c['rejected_at_005'] for c in cases),
                    log_evalue_quantiles=np.quantile([c['log_evalue'] for c in cases],[0,.25,.5,.75,1]).tolist(),
                    max_gap=max(c['fit']['gap'] for c in cases),
                    max_nodes=max(len(c['fit']['nodes']) for c in cases),
                    converged=sum(c['fit']['converged'] for c in cases)))
            record.update(complete=True,scientific_run_complete=True,summary=summary);save()
            summaries.append(dict(dataset=ds,summary=summary,result_sha256=sha(path)))
        except Exception as error:
            record.update(complete=True,scientific_run_complete=False,error=repr(error),elapsed_seconds=time.perf_counter()-started);save();raise
    (out/'summary.json').write_text(json.dumps(dict(complete=True,records=summaries,
        seconds=time.perf_counter()-started,scope='Shared positive scale and unknown finite-catalog viewing law; oracle predictor on reused simulated observations.'),indent=2)+'\n')


if __name__ == '__main__':
    main()
