#!/usr/bin/env python3
"""Summarize both declared prior sweeps without selecting a favorable prior."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fourier_splats.uq_pilot_targets import read_target_lock, LOCK_SHA256

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def main():
    out = BASE/'pilot-selected-fourier-summary-v1'; out.mkdir(exist_ok=True)
    if (out/'summary.json').exists():
        raise RuntimeError('Preserve the previous summary')
    lock = read_target_lock(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
    fig, axes = plt.subplots(2, 3, figsize=(7.15, 4.1), sharex=True, sharey=True)
    source_hashes = {}; records = []; dataset_summaries = []; table_rows = []
    for col, dataset in enumerate(lock['datasets']):
        ds = dataset['dataset']; cases = []
        for directory, expected in [('pilot-selected-fourier-baselines-v1', 12), ('pilot-selected-fourier-baselines-weak', 8)]:
            path = BASE/directory/f'{ds}.json'; raw = path.read_bytes(); result = json.loads(raw)
            if (not result.get('complete') or result.get('error') or len(result['records']) != expected
                    or result['target_lock_sha256'] != LOCK_SHA256 or any(r['numerical_failure'] for r in result['records'])):
                raise ValueError(f'Incomplete or changed baseline family: {path}')
            source_hashes[str(path.relative_to(ROOT))] = hashlib.sha256(raw).hexdigest()
            cases += result['records']
        dataset_summaries.append({'dataset': ds, 'geometry_indices_sha256': result['geometry_indices_sha256'],
            'forward_interpolation_relative_error': result['forward_interpolation_relative_error']})
        priors = sorted({r['prior_deviation_norm'] for r in cases})
        if len(priors) != 5:
            raise ValueError('All five prior scales required')
        for row_index, method in enumerate(['full_fourier_gaussian', 'diagonal_fourier_vi']):
            matrix = np.empty((4, 5))
            for fi, feature in enumerate(dataset['features']):
                for pi, prior in enumerate(priors):
                    case = next(r for r in cases if r['target'] == feature['name'] and r['prior_deviation_norm'] == prior)
                    check = next(c for c in case['checks'] if c['method'] == method and c['scenario'] == 'nominal')
                    matrix[fi, pi] = check['analytic_fixed_signal_coverage']
            ax = axes[row_index, col]
            ax.imshow(matrix, vmin=0, vmax=1, cmap='YlGnBu', aspect='auto')
            for (i, j), value in np.ndenumerate(matrix):
                label = '<.001' if value < .0005 else f'{value:.3f}'
                ax.text(j, i, label, ha='center', va='center', fontsize=6.4, color='white' if value > .65 else 'black')
            ax.axvline(2.5, color='white', linestyle='--', linewidth=1.2)
            ax.set_xticks(range(5), ['0.5', '1', '2', '19.0*', '190*']); ax.tick_params(labelsize=7)
            ax.set_yticks(range(4), ['Region 1', 'Region 2', 'Region 3', 'Center'])
            if row_index == 0:
                ax.set_title(f'EMPIAR-{ds}', fontsize=9)
            else:
                ax.set_xlabel('Prior RMS deviation norm', fontsize=8)
            if col == 0:
                ax.set_ylabel('Full posterior' if row_index == 0 else 'Diagonal VI', fontsize=9)
        for case in cases:
            for check in case['checks']:
                records.append(dict(dataset=ds, feature=case['target'], prior_deviation_norm=case['prior_deviation_norm'],
                    prior_coordinate_sd=case['prior_coordinate_sd'], post_outcome_sensitivity=case['prior_coordinate_sd'] >= .099,
                    **check))
        for prior in priors:
            selected = [r for r in records if r['dataset'] == ds and r['prior_deviation_norm'] == prior
                        and r['scenario'] == 'nominal' and r['method'] == 'full_fourier_gaussian']
            table_rows.append({'dataset': ds, 'prior_deviation_norm': prior,
                'width_range': [min(r['relative_half_width'] for r in selected), max(r['relative_half_width'] for r in selected)],
                'fixed_reference_coverage_range': [min(r['analytic_fixed_signal_coverage'] for r in selected), max(r['analytic_fixed_signal_coverage'] for r in selected)]})
    fig.tight_layout()
    fig.savefig(ROOT/'paper/figures/pilot-selected-fourier.pdf')
    fig.savefig(out/'coverage.png', dpi=200); plt.close(fig)
    with (out/'all-checks.csv').open('w') as file:
        writer = csv.DictWriter(file, fieldnames=records[0], lineterminator='\n'); writer.writeheader(); writer.writerows(records)
    summary = {'complete': True, 'target_lock_sha256': LOCK_SHA256, 'source_hashes': source_hashes,
        'datasets': dataset_summaries, 'nominal_full_posterior_ranges': table_rows,
        'records': records, 'scope': 'Analytic conditional known-noise baseline sensitivity; pointwise reference coverage and prior coverage are different targets.'}
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    tex = [r'\begin{table*}[t]', r'\centering\small', r'\begin{tabular}{lrrr}', r'\toprule',
           r'Stack & Prior RMS norm & Relative half-width & Pointwise coverage\\', r'\midrule']
    for row in table_rows:
        prior = row['prior_deviation_norm']; label = f'{prior:g}' if prior <= 2 else f'{prior:.1f}*'
        widths = row['width_range']; coverage = row['fixed_reference_coverage_range']
        tex.append(f"{row['dataset']} & {label} & {widths[0]:.3f}--{widths[1]:.3f} & {coverage[0]:.3f}--{coverage[1]:.3f}\\\\")
    tex += [r'\bottomrule', r'\end{tabular}',
        r'\caption{Full Fourier-posterior ranges over all four locked targets per stack at nominal poses. Widths use the continuous radius-two no-data scale only for display; the priors are not hard balls. Every prior is retained. Asterisks denote the broader-prior sensitivity specified after the original baseline outcomes. Pointwise coverage integrates prescribed white Gaussian measurement noise for the fixed continuous map generator; it is not experimental-density calibration.}',
        r'\label{tab:pilot-fourier}', r'\end{table*}']
    (ROOT/'paper/tables/pilot-selected-fourier.tex').write_text('\n'.join(tex)+'\n')
    print('COMPLETE', len(records), 'conditional outcome records; 60 matched posterior fits', flush=True)


if __name__ == '__main__':
    main()
