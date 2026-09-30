#!/usr/bin/env python3
"""Summarize the whole declared local-pose baseline without selecting a prior."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'


def main():
    out = BASE/'pilot-selected-fourier-pose-summary-v1'
    if (out/'summary.json').exists():
        raise RuntimeError('Preserve completed summaries')
    source_hashes = {}; rows = []; linearization = []; timings = []; identities = []
    solves = 0; checks = 0; fixed_error = []
    for ds in ['10028', '10049', '10076']:
        path = BASE/f'pilot-selected-fourier-pose-v1/{ds}.json'
        raw = path.read_bytes(); data = json.loads(raw)
        if not data.get('complete') or data.get('error') or len(data['records']) != 16:
            raise ValueError(f'Complete sixteen-case dataset required: {ds}')
        source_hashes[str(path.relative_to(ROOT))] = hashlib.sha256(raw).hexdigest()
        timings.append({'dataset': ds, 'seconds': data['seconds'], 'peak_resident_bytes': data['peak_resident_bytes']})
        linearization += [dict(dataset=ds, **r) for r in data['linearization_checks']]
        for record in data['records']:
            if record['numerical_failure'] or len(record['checks']) != 8:
                raise ValueError('Numerical failure or missing scenario; inspect raw outcomes')
            solves += 1; checks += len(record['checks'])
            identities.append(record['prior_variance_identity_relative_error'])
            if record['fixed_repeat_relative_error'] is not None:
                fixed_error.append(record['fixed_repeat_relative_error'])
        for tau in [.1, 1.]:
            selected = [r for r in data['records'] if r['prior_coordinate_sd'] == tau and r['pose_model'] == 'local_gaussian_pose']
            if len(selected) != 4:
                raise ValueError('Missing target')
            continuous = [r for case in selected for r in case['checks'] if r['scenario'] != 'matched_fourier_model']
            inflation = [100*(r['width_over_original_fixed']-1) for r in selected]
            coverage = [r['analytic_fixed_signal_coverage'] for r in continuous]
            relative_width = [r['relative_half_width'] for r in continuous]
            rows.append({'dataset': ds, 'prior_coordinate_sd': tau,
                'width_increase_percent_min': min(inflation), 'width_increase_percent_max': max(inflation),
                'continuous_reference_coverage_min': min(coverage), 'continuous_reference_coverage_max': max(coverage),
                'continuous_reference_power_min': min(r['correct_sign_probability'] for r in continuous),
                'relative_half_width_min': min(relative_width), 'relative_half_width_max': max(relative_width),
                'matched_prior_predictive_coverage_min': min(r['full_prior_predictive_coverage'] for r in selected),
                'matched_prior_predictive_coverage_max': max(r['full_prior_predictive_coverage'] for r in selected)})
    assert solves == 48 and checks == 384
    summary = {'complete': True, 'source_hashes': source_hashes,
        'scope': 'Local supplied Gaussian pose prior; analytic conditional reference checks, not experimental calibration or Rangan-pipeline reproduction.',
        'solves': solves, 'analytic_reference_checks': checks, 'summary_rows': rows,
        'maximum_variance_identity_relative_error': max(identities),
        'maximum_fixed_control_repeat_relative_error': max(fixed_error),
        'linearization_checks': linearization, 'runtime': timings,
        'interpretation': 'Small width inflation describes the fixed pilot Jacobian model only. Gaussian coordinate SDs are not hard joint radii; density-pose products and nonlinear terms are omitted.'}
    out.mkdir(exist_ok=True)
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    table = [r'\begin{table*}[t]', r'\centering',
        r'\caption{Local Gaussian pose baseline on the twelve locked features. Width increase is relative to the matching fixed-pose full posterior. Coverage ranges and minimum sign power use the continuous reference at nominal and the prescribed boundary poses, with known measurement noise. They are pointwise sensitivity results, not empirical density calibration. Nominal marginal prior coverage is $0.995833$.}',
        r'\label{tab:fourier-pose}', r'\begin{tabular}{llrrr}', r'\toprule',
        r'Stack & Prior coordinate SD & Width increase (\%) & Reference coverage & Min. sign power\\', r'\midrule']
    for r in rows:
        table.append(f"{r['dataset']} & {r['prior_coordinate_sd']:g} & {r['width_increase_percent_min']:.4f}--{r['width_increase_percent_max']:.4f} & {r['continuous_reference_coverage_min']:.4f}--{r['continuous_reference_coverage_max']:.4f} & {r['continuous_reference_power_min']:.4f}\\\\")
    table += [r'\bottomrule', r'\end{tabular}', r'\end{table*}']
    (ROOT/'paper/tables/pilot-fourier-pose.tex').write_text('\n'.join(table)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
