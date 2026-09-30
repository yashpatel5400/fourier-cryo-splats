#!/usr/bin/env python3
"""Summarize all existing matched-baseline pose checks without new fitting."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    source = ROOT/'results/uncertainty/development/pilot-selected-fourier-summary-v1/summary.json'
    raw = source.read_bytes(); old = json.loads(raw)
    if not old.get('complete') or len(old['records']) != 960:
        raise ValueError('Require the entire completed five-prior baseline family')
    out = source.parent/'pose-sensitivity-summary.json'
    if out.exists():
        raise RuntimeError('Preserve the existing derived summary')
    for name, expected in old['source_hashes'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Original baseline evidence changed')
    records = old['records']
    priors = sorted({r['prior_deviation_norm'] for r in records})
    methods = sorted({r['method'] for r in records})
    scenarios = sorted({r['scenario'] for r in records})
    result = []
    for prior in priors:
        for method in methods:
            for scenario in scenarios:
                rows = [r for r in records if r['prior_deviation_norm'] == prior
                        and r['method'] == method and r['scenario'] == scenario]
                if len(rows) != 12 or len({(r['dataset'], r['feature']) for r in rows}) != 12:
                    raise ValueError('Missing or duplicate target outcomes')
                result.append(dict(prior_deviation_norm=prior,
                    prior_coordinate_sd=rows[0]['prior_coordinate_sd'],
                    post_outcome_sensitivity=rows[0]['post_outcome_sensitivity'],
                    method=method, scenario=scenario, targets=12,
                    minimum_coverage=min(r['analytic_fixed_signal_coverage'] for r in rows),
                    maximum_coverage=max(r['analytic_fixed_signal_coverage'] for r in rows),
                    minimum_sign_power=min(r['correct_sign_probability'] for r in rows),
                    maximum_sign_power=max(r['correct_sign_probability'] for r in rows),
                    minimum_relative_half_width=min(r['relative_half_width'] for r in rows),
                    maximum_relative_half_width=max(r['relative_half_width'] for r in rows)))
    if len(result) != 80:
        raise ValueError('Expected all five priors, two methods and eight scenarios')
    outcome = {'complete': True, 'scope': 'Descriptive regrouping of existing outcomes, not additional experiments or experimental coverage',
        'source': str(source.relative_to(ROOT)), 'source_sha256': hashlib.sha256(raw).hexdigest(),
        'runner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'records': result, 'nominal_marginal_level': 1-.05/12,
        'caution': 'These are specified map-like generators and two pose patterns, not a worst-case pose guarantee. Fixed-signal coverage differs from prior-predictive coverage. The broader two priors remain post-outcome sensitivities.'}
    out.write_text(json.dumps(outcome, indent=2)+'\n')
    with out.with_suffix('.csv').open('w') as file:
        writer = csv.DictWriter(file, fieldnames=result[0], lineterminator='\n')
        writer.writeheader(); writer.writerows(result)
    print('COMPLETE', len(records), 'existing checks grouped into', len(result), 'descriptive rows', flush=True)


if __name__ == '__main__':
    main()
