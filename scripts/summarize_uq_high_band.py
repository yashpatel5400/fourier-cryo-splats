#!/usr/bin/env python3
"""Retain all completed higher-band continuous targets, including low sign power."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
records = []; sources = {}
for dataset in ['10028', '10049', '10076']:
    folder = 'continuous-quadrature-high-band-probe' if dataset == '10028' else 'continuous-quadrature-high-band-additional'
    path = BASE/folder/f'{dataset}.json'; data = json.loads(path.read_text())
    sources[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    if len(data['targets']) != 2 or data['config']['frequency_radius'] != 12:
        raise AssertionError('Higher-band study incomplete')
    for row in data['targets']:
        fit = row['fit']; fine = next(r for r in row['cell_projection_diagnostics'] if r['box'] == 64)
        records.append({'dataset': dataset, 'target': row['target'], 'relative_half_width': row['width_fraction_of_no_data'],
                        'relative_gap': max(0., (fit['objective']-fit['dual_lower_bound'])/fit['objective']),
                        'converged': fit['converged'], 'solve_and_audit_seconds': row['seconds'],
                        'setup_seconds_shared_across_targets': data['setup_seconds'],
                        'stored_arrays_bytes': data['stored_gram_or_quadrature_bytes'],
                        'reference_correct_sign_probability': row['reference_correct_sign_probability'],
                        'cell64_width_over_continuous_width': fine['width_fraction_of_continuous'],
                        'cell64_coverage_at_continuous_boundary': fine['coverage_at_continuous_bias_boundary'],
                        'reference_coverage': row['reference_analytic_coverage']})
summary = {'stage': 'development higher-band continuous fixed-pose audit', 'sources': sources,
           'analysis_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'settings': len(records), 'records': records,
           'interpretation': ['128 particles; Fourier radius 12; Gaussian feature width 0.03 of field; fixed poses.',
                              'A narrower interval than no-data does not imply power for a particular biological feature.',
                              'Continuous quadrature error is bounded analytically; dense independent sinc check is skipped at this scale.',
                              'Timings share concurrent Mac load and are not isolated performance measurements.']}
out = BASE/'continuous-high-band-summary'; out.mkdir(exist_ok=True)
(out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
