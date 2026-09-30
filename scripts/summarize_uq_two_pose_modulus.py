#!/usr/bin/env python3
"""Summarize the complete prespecified two-pose grid, retaining failed attempts."""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    folder=BASE/'two-pose-modulus-v3';rows=[];sources={};attempts=[]
    cases=[('nominal',0),('coherent_x',1),('random_boundary',1),('coherent_x',2),('random_boundary',2)]
    for dataset in ['10028','10049','10076']:
        for target in ['center','contrast']:
            for scenario,degrees in cases:
                path=folder/f'{dataset}-{target}-{scenario}-{degrees}.json';d=json.loads(path.read_text())
                if not d.get('complete') or d.get('error') or d.get('numerical_failure'):raise ValueError('Require the entire successful v3 grid')
                digest=hashlib.sha256(path.with_suffix('.npz').read_bytes()).hexdigest()
                assert digest==d['witness_array_sha256']
                sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
                fit=d['fit'];witness=fit['witness'];normalizer=d['density_radius']*fit['target_norm']
                assert abs(normalizer-d['no_data_half_width'])<1e-12
                assert d['maximum_scaled_pose_norm']<=1+1e-12 and d['minimum_sum_testing_error']>2*d['alpha']
                rows.append({'dataset':dataset,'target':target,'scenario':scenario,'degrees':degrees,
                    'relative_half_width_lower':d['relative_half_width_lower'],
                    'half_width_lower':witness['fixed_length_half_width_lower'],'no_data_half_width':normalizer,
                    'fixed_pair_modulus_gap':fit['relative_modulus_gap'],'converged':fit['converged'],
                    'minimum_sum_testing_error':d['minimum_sum_testing_error'],'seconds':d['seconds'],
                    'best_lower_iteration':fit['best_lower_iteration'],'best_upper_iteration':fit['best_upper_iteration'],
                    'source':str(path.relative_to(ROOT))})
    for directory in sorted(BASE.glob('two-pose-modulus*')):
        for path in sorted(directory.glob('*.json')):
            if path.name.endswith('-setup.json'):continue
            d=json.loads(path.read_text())
            if not d.get('stage','').startswith('Constructive conditional two-pose'):continue
            attempts.append({'source':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                'complete':bool(d.get('complete') and not d.get('error')),'error':d.get('error')})
    result={'complete':True,'scope':'Selected-pose lower bounds for deterministic-length intervals under the declared continuous L2 ball and white Gaussian noise. Not a global pose modulus or experimental calibration.',
        'primary_grid':'two-pose-modulus-v3','source_hashes':sources,'records':rows,'all_variant_attempts':attempts,
        'primary_cases':len(rows),'maximum_fixed_pair_relative_gap':max(r['fixed_pair_modulus_gap'] for r in rows),
        'translation_radius_A':.5,'density_radius':2.,'target_width_fraction_field':.07,
        'pose_set':'Per-particle joint five-dimensional unit ball with antipodal poses',
        'normalization':'Lower half-width divided by B||ell||, matching the upper-audit no-data half-width; unnormalized values are retained.'}
    out=BASE/'two-pose-modulus-summary';out.mkdir(exist_ok=True)
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    with (out/'cases.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n');writer.writeheader();writer.writerows(rows)
    lines=[r'\begin{table*}[t]',r'\centering',r'\begin{tabular}{llrrrrr}',r'\toprule',
           r'Stack & Target & Nominal & Coherent $1^\circ$ & Random $1^\circ$ & Coherent $2^\circ$ & Random $2^\circ$ \\',r'\midrule']
    for dataset in ['10028','10049','10076']:
        for target in ['center','contrast']:
            values=[next(r['relative_half_width_lower'] for r in rows if r['dataset']==dataset and r['target']==target and r['scenario']==scenario and r['degrees']==degrees) for scenario,degrees in cases]
            lines.append(f'{dataset} & {target} & '+' & '.join(f'{v:.3f}' for v in values)+r' \\')
    lines.extend([r'\bottomrule',r'\end{tabular}',
        r'\caption{Constructive lower bounds on the half-width of any uniformly valid deterministic-length interval, divided by the no-data half-width $B\norm\ell$. Each entry uses two feasible antipodal pose configurations, continuous densities in the original ball and a Gaussian testing distance below $2\Phi^{-1}(0.95)$. Random poses are fixed design draws. All 30 fixed-pair modulus gaps are below 0.5\%; this does not certify maximization over poses. A small entry cannot establish that a larger upper bound is loose.}',
        r'\label{tab:two-pose-modulus}',r'\end{table*}'])
    (ROOT/'paper/tables/two-pose-modulus.tex').write_text('\n'.join(lines)+'\n')
    print('DONE',len(rows),result['maximum_fixed_pair_relative_gap'])


if __name__=='__main__':main()
