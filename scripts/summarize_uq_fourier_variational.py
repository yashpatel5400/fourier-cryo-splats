#!/usr/bin/env python3
"""Create transparent baseline comparisons from all completed prior sweeps."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def main():
    out = BASE/'fourier-variational-summary'; out.mkdir(exist_ok=True)
    figure, axes = plt.subplots(2, 3, figsize=(7.1, 4.25), sharex=True, sharey=True)
    hashes = {}; records = []; datasets = []
    tex = [r'\begin{table*}[t]', r'\centering\small', r'\begin{tabular}{llrrr}', r'\toprule',
           r'Stack & Target & VI/full width & Full coverage & VI coverage\\', r'\midrule']
    for column, dataset in enumerate(['10028', '10049', '10076']):
        path = BASE/'fourier-variational-baseline'/f'{dataset}.json'
        d = json.loads(path.read_text())
        if not d['complete'] or len(d['records']) != 6:
            raise RuntimeError('Incomplete baseline; do not omit failed priors')
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        datasets.append({k: d[k] for k in ['dataset', 'particles', 'parameters', 'seconds',
            'sparse_matrix_bytes', 'peak_resident_bytes', 'continuous_reference_relative_forward_error']})
        for row, target in enumerate(['center', 'contrast']):
            ax = axes[row, column]; cases = [r for r in d['records'] if r['target'] == target]
            for method, color, marker, label in [('full_fourier_gaussian', '#246a9b', 'o', 'Full posterior'),
                                                ('diagonal_fourier_vi', '#b54c36', 's', 'Diagonal VI')]:
                coverage = []
                for case in cases:
                    check = next(c for c in case['checks'] if c['method'] == method and c['scenario'] == 'continuous_reference_nominal')
                    coverage.append(check['analytic_fixed_signal_coverage'])
                ax.plot([np.sqrt(c['prior_expected_squared_deviation_norm']) for c in cases], coverage,
                        color=color, marker=marker, label=label)
            ax.axhline(.95, color='#777777', linewidth=.7, linestyle='--')
            ax.set_xscale('log', base=2); ax.set_xticks([.5, 1, 2], ['0.5', '1', '2'])
            ax.set_ylim(-.025, 1.025); ax.grid(alpha=.15); ax.tick_params(labelsize=8)
            if row == 0: ax.set_title('EMPIAR-'+dataset, fontsize=9)
            if row == 1: ax.set_xlabel('Prior deviation RMS norm', fontsize=8)
            if column == 0: ax.set_ylabel(target.capitalize()+' coverage', fontsize=8)
            if row == 0 and column == 2: ax.legend(fontsize=7, loc='center right', frameon=False)
            for case in cases:
                checks = {c['method']: c for c in case['checks'] if c['scenario'] == 'continuous_reference_nominal'}
                records.append({'dataset': dataset, 'target': target,
                    'prior_expected_squared_deviation_norm': case['prior_expected_squared_deviation_norm'],
                    'width_ratio': case['vi_over_full_posterior_width'],
                    'full_prior_coverage': case['full_prior_predictive_coverage'],
                    'vi_prior_coverage': case['diagonal_vi_prior_predictive_coverage'],
                    'identity_relative_error': case['prior_predictive_variance_identity_relative_error'],
                    'continuous_nominal_checks': checks})
                if case['prior_expected_squared_deviation_norm'] == 1:
                    full = checks['full_fourier_gaussian']['analytic_fixed_signal_coverage']
                    vi = checks['diagonal_fourier_vi']['analytic_fixed_signal_coverage']
                    tex.append(f"{dataset} & {target.capitalize()} & {case['vi_over_full_posterior_width']:.3f} & {full:.3f} & {vi:.3f}\\\\")
    figure.tight_layout()
    figure.savefig(ROOT/'paper/figures/fourier-variational.pdf')
    figure.savefig(out/'coverage.png', dpi=170); plt.close(figure)
    tex += [r'\bottomrule', r'\end{tabular}',
        r'\caption{Fixed-pose Fourier Gaussian baselines at prior expected squared deviation norm one, 1,024 particles per stack, radius 12 and 10\,\AA\ target sigma. Coverage integrates known measurement noise for the continuous deposited-map generator at nominal poses; it is not empirical experimental coverage. The full posterior has exactly 0.95 coverage under its matched prior predictive distribution. VI is wider for these targets. Poor fixed-map coverage includes shrinkage and interpolation/representation error; it does not refute the source method\textquotesingle s stated probabilistic assumptions.}',
        r'\label{tab:fouriervi}', r'\end{table*}']
    (ROOT/'paper/tables/fourier-variational.tex').write_text('\n'.join(tex)+'\n')
    result = {'complete': True, 'source_sha256': hashes, 'datasets': datasets, 'records': records,
              'width_ratio_range': [min(r['width_ratio'] for r in records), max(r['width_ratio'] for r in records)],
              'vi_prior_coverage_range': [min(r['vi_prior_coverage'] for r in records), max(r['vi_prior_coverage'] for r in records)],
              'maximum_variance_identity_relative_error': max(r['identity_relative_error'] for r in records)}
    audits = []
    for dataset in ['10028', '10049', '10076']:
        path = BASE/'fourier-variational-continuous-audit'/f'{dataset}.json'
        d = json.loads(path.read_text())
        if not d['complete'] or len(d['records']) != 6:
            raise RuntimeError('Incomplete continuous re-audit')
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        audits.extend(d['records'])
    result['continuous_same_estimator_audit'] = {
        'cases': len(audits),
        'relative_width_range': [min(r['relative_half_width_before_fallback'] for r in audits),
                                 max(r['relative_half_width_before_fallback'] for r in audits)],
        'posterior_width_inflation_range': [min(r['half_width_over_full_posterior'] for r in audits),
                                           max(r['half_width_over_full_posterior'] for r in audits)],
        'minimum_reference_coverage': min(r['reference_same_estimator_coverage'] for r in audits),
        'maximum_reference_power': max(r['reference_same_estimator_correct_sign_probability'] for r in audits),
        'maximum_recentered_reference_power': max(r['reference_recentered_correct_sign_probability'] for r in audits),
        'maximum_reference_bias_fraction': max(abs(r['reference_expected_center']-r['reference_target'])/r['same_estimator_bias_upper'] for r in audits)}
    (out/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['width_ratio_range', 'vi_prior_coverage_range', 'maximum_variance_identity_relative_error']}, indent=2))


if __name__ == '__main__': main()
