#!/usr/bin/env python3
"""Post-audit bounds, frequency checks and amplitude-mass sensitivity."""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fourier_splats.uq_covariance_cone_bounds import (scalar_dual_tangent_upper,
    unrestricted_matrix_growth_upper)

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    out=BASE/'paired-statistics-bound-corrections-v1'
    if out.exists():raise RuntimeError('Preserve previous correction')
    out.mkdir();record=dict(power_cases=[],matrix_cases=[],mass_sensitivity=[],input_hashes={})
    def source(path):record['input_hashes'][str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    p=BASE/'paired-power-enlarged-cone-v1/summary.json';power=json.loads(p.read_text());source(p)
    p=BASE/'paired-covariance-finite-view-v1/summary.json';matrix=json.loads(p.read_text());source(p)
    for ds in ['10028','10049','10076']:
        p=ROOT/power['arrays'][ds]['path'];source(p)
        with np.load(p) as d:
            q=d['q'];keys={tuple(v) for v in q}
            assert len(keys)==len(q) and (0.,0.) not in keys and not any(tuple(-v) in keys for v in q) and abs(q).max()<32
            signal=d['signal']
            for case in power['cases']:
                if case['dataset']!=ds:continue
                key=f"{case['candidate']}_{case['added_views']}"
                approximation=d[key+'_approximation'];values={}
                for v in [1.,2.,4.]:
                    bound=sum(scalar_dual_tangent_upper(s*tr/v,a*tr/v)['upper']
                        for transfer in d['transfer_squared'] for s,a,tr in zip(signal,approximation,transfer))*16
                    values[str(v)]=bound
                record['power_cases'].append(dict(dataset=ds,candidate=case['candidate'],views=case['total_views'],
                    corrected_expected_log_128_upper=values,nonredundant_frequencies=True))
            if ds=='10076':
                for candidate,name in [(1,'region_half_removed'),(2,'region_removed')]:
                    powers=abs(d['fourier'][candidate])**2;scale=np.maximum(signal,1e-6*signal.max())
                    a=powers.T/scale[:,None];b=signal/scale
                    for sign,label in [(1.,'minimum'),(-1.,'maximum')]:
                        result=linprog(sign*np.ones(len(powers)),A_eq=a,b_eq=b,bounds=(0,None),method='highs',
                            options=dict(time_limit=60.,primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
                        row=dict(dataset=ds,candidate=name,objective=label,status=int(result.status),message=result.message)
                        if result.success:
                            weights=np.maximum(result.x,0.);np.savez_compressed(out/f'{ds}-{name}-{label}-mass.npz',weights=weights)
                            row.update(mass=float(weights.sum()),maximum_scaled_residual=float(np.max(abs(a@weights-b))))
                        record['mass_sensitivity'].append(row)
        p=ROOT/matrix['arrays'][ds]['path'];source(p)
        with np.load(p) as d:
            for case in matrix['cases']:
                if case['dataset']!=ds:continue
                key=f"{case['candidate']}_{case['rank']}"
                delta=d[key+'_signal_second_moment']-d[key+'_approximation']
                bounds={str(v):unrestricted_matrix_growth_upper(delta,v) for v in [1.,2.,4.]}
                record['matrix_cases'].append(dict(dataset=ds,candidate=case['candidate'],rank=case['rank'],
                    spectral_growth_uppers=bounds,
                    lp_primal_dual_gap=abs(case['lp_objective']-case['dual_objective']),
                    original_repair=case['direction_repair']))
    record['complete']=True
    record['scope']='Root-error-corrected numerical bounds; no validated arithmetic, continuous certification, or new observations.'
    (out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(power_cases=len(record['power_cases']),matrix_cases=len(record['matrix_cases']),mass=record['mass_sensitivity']),indent=2))


if __name__=='__main__':main()
