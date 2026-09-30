#!/usr/bin/env python3
"""Local cell-quality diagnostic, not a global reconstruction/validation run."""
from itertools import product
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from scipy.optimize import minimize
from fourier_splats.uq_mixture_curvature import CurvatureGaussianOrbit
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
PROTOCOL=ROOT/'research/uncertainty/CONTINUOUS-MIXTURE-CURVATURE-PROTOCOL.md'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    files=[Path(__file__),PROTOCOL,
        ROOT/'research/uncertainty/CONTINUOUS-MIXTURE-CURVATURE-THEORY.md',
        ROOT/'src/fourier_splats/uq_mixture_curvature.py',
        ROOT/'src/fourier_splats/uq_continuous_mixture.py',
        ROOT/'tests/test_uq_mixture_curvature.py']
    for p in files:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise ValueError('Commit source/protocol before case computation')
    out=BASE/'continuous-mixture-curvature-v1'
    if out.exists():raise RuntimeError('Preserve every previous outcome')
    out.mkdir();summaries=[];calculation_seconds=0.
    for ds in ['10028','10049','10076']:
        path=out/f'{ds}.json'
        record=dict(complete=False,dataset=ds,cases=[],
            source_snapshot=source_snapshot(ROOT,Path(__file__),[str(p.relative_to(ROOT)) for p in files[1:]]),
            scope='Local bound diagnosis at oracle simulated orientations; not a global SO(3) calculation or experimental validation.')
        def save():path.write_text(json.dumps(record,indent=2)+'\n')
        save()
        try:
            source=BASE/f'continuous-mixture-v1/{ds}.json';deadline=time.monotonic()+3600
            while True:
                try:ready=source.exists() and json.loads(source.read_text()).get('complete')
                except json.JSONDecodeError:ready=False
                if ready:break
                if time.monotonic()>deadline:raise TimeoutError('Parent wait budget reached')
                time.sleep(10)
            parent=json.loads(source.read_text())
            if not parent.get('scientific_run_complete'):raise ValueError('Completed parent required')
            array_path=ROOT/parent['arrays_file']
            if sha(array_path)!=parent['arrays_sha256']:raise ValueError('Parent arrays changed')
            arrays=np.load(array_path)
            record['input_hashes']={str(source.relative_to(ROOT)):sha(source),str(array_path.relative_to(ROOT)):sha(array_path)}
            for position in [0,64,127]:
                orbit=CurvatureGaussianOrbit(arrays['plane'],arrays['centers'],arrays['sigma'],
                    arrays['coefficients'],arrays['transfer'][position:position+1],arrays['observed'][position:position+1])
                center=arrays['oracle_angles'][position]
                for budget in [.1,.5,1.,2.,5.,10.]:
                    if calculation_seconds>1800:raise TimeoutError('Thirty-minute calculation budget reached')
                    begin=time.perf_counter();half=np.full(3,np.deg2rad(budget)/3)
                    lower,upper=center-half,center+half
                    cell=orbit.cell(lower,upper)
                    def objective(angle):
                        mean,jac,_=orbit.mean_jacobian(angle);residual=orbit.observed-mean
                        return .5*float(np.sum(residual**2)),-np.einsum('nd,ndj->j',residual,jac)
                    starts=[center]+[center+half*np.array(s) for s in product([-1,1],repeat=3)]
                    local=[]
                    for initial in starts:
                        result=minimize(objective,initial,jac=True,method='L-BFGS-B',
                            bounds=list(zip(lower,upper)),options=dict(maxiter=100,ftol=1e-12,gtol=1e-9))
                        local.append(dict(initial=initial.tolist(),angles=result.x.tolist(),
                            log_kernel=-float(result.fun),success=bool(result.success),status=int(result.status),
                            message=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev)))
                    feasible=max(float(cell['exact'][0]),max(v['log_kernel'] for v in local))
                    bound=float(cell['envelope'][0])
                    case=dict(position=position,budget_degrees=budget,center=center.tolist(),half_widths_radians=half.tolist(),
                        center_log_kernel=float(cell['exact'][0]),original_upper=float(cell['first_order_envelope'][0]),
                        quadratic_upper=float(cell['quadratic_envelope'][0]),combined_upper=bound,
                        third_log_remainder=float(cell['third_log_remainder'][0]),
                        quadratic_increment_upper=float(cell['quadratic_increment_upper'][0]),
                        feasible_log_kernel=feasible,original_gap=float(cell['first_order_envelope'][0])-feasible,
                        combined_gap=bound-feasible,local_searches=local,seconds=time.perf_counter()-begin)
                    record['cases'].append(case);save();print(ds,position,budget,case['original_gap'],case['combined_gap'],flush=True)
                    if bound<feasible-1e-7:raise ArithmeticError('Cell bound violated by feasible search')
                    calculation_seconds+=case['seconds']
            record.update(complete=True,scientific_run_complete=True);save()
            summary=[]
            for budget in [.1,.5,1.,2.,5.,10.]:
                cases=[c for c in record['cases'] if c['budget_degrees']==budget]
                summary.append(dict(budget_degrees=budget,cases=len(cases),
                    original_gap_median=float(np.median([c['original_gap'] for c in cases])),
                    combined_gap_median=float(np.median([c['combined_gap'] for c in cases])),
                    combined_gap_max=max(c['combined_gap'] for c in cases)))
            summaries.append(dict(dataset=ds,summary=summary,result_sha256=sha(path)))
        except Exception as error:
            record.update(complete=True,scientific_run_complete=False,error=repr(error));save();raise
    (out/'summary.json').write_text(json.dumps(dict(complete=True,records=summaries,
        calculation_seconds=calculation_seconds,scope='Fifty-four oracle-centered local envelope checks; no global optimality or density calibration.'),indent=2)+'\n')


if __name__=='__main__':main()
