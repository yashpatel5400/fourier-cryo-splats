#!/usr/bin/env python3
"""Summarize every case in the locked continuous-v1 study, only when complete.

This analysis does not modify the frozen estimator, source lock, or outcomes.
Cases sharing a geometry or paired noise draws are not treated as independent
biological replications. Exact coverage is the primary noise-marginal quantity.
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--study', choices=['continuous-v1', 'continuous-moments-v2'], default='continuous-v1')
ARGS = parser.parse_args()
BASE = ROOT/'results/uncertainty/confirmation'/ARGS.study
OUT = BASE/'summary'
LOCK = ROOT/'research/uncertainty/confirmation'/ARGS.study/'locked-study.json'
METHODS = ['noise_only', 'cell24_fixed_pose', 'continuous_fixed_pose',
           'continuous_pose', 'no_data']
LABELS = ['Noise only', '24-cell audit', 'Continuous, fixed pose',
          'Continuous, bounded pose', 'No data']
if ARGS.study == 'continuous-moments-v2':
    METHODS.insert(3, 'continuous_pose_spatial_maximum')
    LABELS.insert(3, 'Pose, spatial maximum')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def range_summary(values):
    x = np.asarray(values, dtype=float)
    return {'count': len(x), 'minimum': float(x.min()), 'median': float(np.median(x)),
            'maximum': float(x.max())}


def main():
    lock = json.loads(LOCK.read_text()); lock_hash = sha(LOCK)
    cfg = lock['config']
    settings_per_geometry = len(cfg['targets'])*len(cfg['widths'])*len(cfg['angles'])
    geometry_count = len(lock['datasets'])*lock['replicates']
    setting_count = settings_per_geometry*geometry_count
    case_count = 24*setting_count
    fit_count = len(cfg['targets'])*len(cfg['widths'])*geometry_count
    sources = []
    for dataset in lock['datasets']:
        for replicate in range(lock['replicates']):
            path = BASE/f'{dataset}-replicate{replicate}.json'
            if not path.exists():
                raise RuntimeError(f'Study incomplete: missing {path.name}')
            data = json.loads(path.read_text())
            if not data['complete'] or data['protocol_lock_sha256'] != lock_hash:
                raise RuntimeError(f'Incomplete or mismatched source: {path.name}')
            if len(data['records']) != settings_per_geometry or any(len(r['cases']) != 24 for r in data['records']):
                raise AssertionError('Unexpected locked-study record count')
            sources.append((path, data))
    OUT.mkdir(exist_ok=True)
    cases = []; designs = []; fits = []; runtimes = []
    for path, data in sources:
        runtimes.append(data['seconds'])
        seen = set()
        for record in data['records']:
            common = {k: record[k] for k in ['target', 'width', 'angle_degrees']}
            common.update(dataset=data['dataset'], replicate=data['replicate'])
            key = record['target'], record['width']
            if key not in seen:
                seen.add(key)
                fit = record['fit']
                fits.append(dict(common, converged=fit['converged'],
                                 relative_gap=max(0., (fit['objective']-fit['dual_lower_bound'])/fit['objective']),
                                 outer_iterations=len(fit['history'])))
            for method in METHODS:
                half = record['half_widths_before_fallback'][method]
                no_data = record['no_data_half_width']
                designs.append(dict(common, method=method,
                                    relative_half_width=min(half/no_data, 1.),
                                    uses_no_data=method == 'no_data' or half >= no_data,
                                    field_gram_bytes=record['audit']['stored_field_and_gram_bytes']))
            for case in record['cases']:
                group = 'continuous_boundary' if case['generator'].startswith('continuous_') else 'cell_generator'
                for method, result in case['methods'].items():
                    cases.append(dict(common, method=method, generator=case['generator'],
                                      generator_group=group, scenario=case['scenario'],
                                      within_declared_model=case['within_declared_model'],
                                      true_target=case['true_target'],
                                      analytic_coverage=result['analytic_coverage'],
                                      monte_carlo_coverage=result['monte_carlo_coverage'],
                                      binomial_low=result['binomial_95'][0],
                                      binomial_high=result['binomial_95'][1],
                                      noise_replicates=result['replications'],
                                      half_width=result['half_width'],
                                      relative_half_width=result['half_width']/no_data,
                                      uses_no_data=result['uses_no_data'],
                                      actual_bias=result['actual_bias'],
                                      actual_noise_sd=result['actual_noise_sd'],
                                      correct_sign_probability=result['correct_sign_probability']))
    if len(cases) != case_count*len(METHODS) or len(designs) != setting_count*len(METHODS) or len(fits) != fit_count:
        raise AssertionError('Do not summarize an incomplete study')
    for name, rows in [('case-results.csv', cases), ('design-results.csv', designs), ('fit-results.csv', fits)]:
        with (OUT/name).open('w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    strata = defaultdict(list)
    for row in cases:
        strata[(row['method'], row['within_declared_model'], row['generator_group'], row['scenario'])].append(row)
    coverage = []
    for (method, inside, group, scenario), rows in strata.items():
        exact = [r['analytic_coverage'] for r in rows]
        coverage.append({'method': method, 'within_declared_model': inside,
                         'generator_group': group, 'scenario': scenario,
                         'coverage': range_summary(exact),
                         'under_95_count': sum(x < .95-1e-8 for x in exact),
                         'sign_probability': range_summary([r['correct_sign_probability'] for r in rows]),
                         'maximum_absolute_mc_error': max(abs(r['monte_carlo_coverage']-r['analytic_coverage']) for r in rows),
                         'analytic_outside_marginal_binomial95_count': sum(not r['binomial_low'] <= r['analytic_coverage'] <= r['binomial_high'] for r in rows)})
    width_groups = defaultdict(list)
    for row in designs:
        width_groups[(row['dataset'], row['target'], row['width'], row['angle_degrees'], row['method'])].append(row)
    widths = []
    for (dataset, target, width, angle, method), rows in width_groups.items():
        widths.append({'dataset': dataset, 'target': target, 'width': width,
                       'angle_degrees': angle, 'method': method,
                       'relative_half_width': range_summary([r['relative_half_width'] for r in rows]),
                       'no_data_design_count': sum(r['uses_no_data'] for r in rows)})
    ratios = []
    for _, data in sources:
        for record in data['records']:
            boundaries = [c for c in record['cases'] if c['generator'].startswith('continuous_pose_boundary')]
            # Noise-only never falls back in this locked design: its actual_bias
            # belongs to the original linear estimator, even if another method
            # would switch to the no-data center.
            if any(c['methods']['noise_only']['uses_no_data'] for c in boundaries):
                raise AssertionError('Need original estimator bias for feasible lower bound')
            ratios.append(max(abs(c['methods']['noise_only']['actual_bias']) for c in boundaries)/record['audit']['total_bias'])
    full_rows = [r for r in cases if r['method'] == 'continuous_pose' and r['within_declared_model']]
    if len(full_rows) != 16*setting_count or min(r['analytic_coverage'] for r in full_rows) < .95-1e-8:
        raise AssertionError('Primary correctness criterion failed; retain and investigate raw outcomes')
    summary = {'status': 'complete frozen-protocol simulation study; conditional acquisition geometries',
               'study': ARGS.study,
               'protocol_lock_sha256': lock_hash,
               'analysis_source_sha256': sha(Path(__file__)),
               'source_files': {str(p.relative_to(ROOT)): sha(p) for p, _ in sources},
               'audit_settings': setting_count, 'signal_scenario_records': case_count,
               'in_class_records': 16*setting_count, 'out_of_class_records': 8*setting_count,
               'unique_fixed_pose_fits': fit_count,
               'minimum_in_class_continuous_pose_coverage': min(r['analytic_coverage'] for r in full_rows),
               'fit_converged_count': sum(r['converged'] for r in fits),
               'fit_relative_gaps': range_summary([r['relative_gap'] for r in fits]),
               'geometry_run_seconds': range_summary(runtimes),
               'total_geometry_run_seconds': sum(runtimes),
               'stored_pose_field_gram_bytes': range_summary([r['field_gram_bytes'] for r in designs]),
               'coherent_pose_feasible_bias_over_uniform_bound': range_summary(ratios),
               'coverage_strata': coverage, 'width_strata': widths,
               'interpretation': ['Analytic coverage integrates known Gaussian measurement noise conditional on each fixed signal/design.',
                                  'Monte Carlo draws are paired across methods per case; marginal binomial checks are not simultaneous intervals.',
                                  'Geometry replicates can share source exposure groups; they are not independent biological datasets.',
                                  'Density radius, support, pose limits and noise scale are supplied simulation assumptions.',
                                  'Feasible bias ratios use one coherent pose, not a global nonlinear optimization.',
                                  'Stored-array memory excludes library temporaries and is not peak resident memory.']}
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    plot(designs, cases, cfg)
    table(widths, setting_count, cfg)
    print(json.dumps({k: summary[k] for k in ['audit_settings', 'signal_scenario_records',
          'minimum_in_class_continuous_pose_coverage', 'fit_converged_count',
          'total_geometry_run_seconds', 'coherent_pose_feasible_bias_over_uniform_bound']}, indent=2))


def plot(designs, cases, cfg):
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 3, figsize=(10., 3.1))
    panels = [(width, None) for width in cfg['widths']] if len(cfg['widths']) == 2 else [(cfg['widths'][0], name) for name in cfg['targets']]
    for j, (width, selected_target) in enumerate(panels):
        ax = axes[j]
        for dataset, color in zip(['10028', '10049', '10076'], ['#2166ac', '#d95f02', '#1b9e77']):
            for target, marker, shift in [('center', 'o', -.025), ('contrast', '^', .025)]:
                if selected_target is not None and target != selected_target:
                    continue
                for i, angle in enumerate(cfg['angles']):
                    values = [r['relative_half_width'] for r in designs if r['method'] == 'continuous_pose'
                              and r['dataset'] == dataset and r['target'] == target
                              and r['width'] == width and r['angle_degrees'] == angle]
                    offset = ['10028', '10049', '10076'].index(dataset)*.14-.14+shift
                    ax.plot(np.full(len(values), i+offset), values, marker, color=color, markersize=3, alpha=.8)
        ax.axhline(1, color='.4', linestyle=':', linewidth=.8)
        ax.set_xticks([0, 1], [f'{angle:g}' for angle in cfg['angles']]); ax.set_xlim(-.4, 1.4); ax.set_ylim(0, 1.04)
        title = f'Target scale {width:.2f} of field' if selected_target is None else ('Central average' if selected_target == 'center' else 'Axial contrast')
        ax.set_title(title); ax.set_xlabel('Rotation radius (degrees)')
    axes[0].set_ylabel('Half-width / no-data half-width')
    groups = [('cell_generator', 'nominal'), ('cell_generator', 'random_boundary'),
              ('cell_generator', 'coherent_boundary'), ('continuous_boundary', 'nominal'),
              ('continuous_boundary', 'coherent_boundary')]
    matrix = np.array([[min(r['analytic_coverage'] for r in cases if r['method'] == method
                           and r['generator_group'] == group and r['scenario'] == scenario
                           and r['within_declared_model']) for group, scenario in groups] for method in METHODS])
    axes[2].imshow(matrix, vmin=0, vmax=1, cmap='viridis', aspect='auto')
    short_labels = ['Noise', 'Cell', 'Fixed', 'Pose', 'No data'] if len(METHODS) == 5 else ['Noise', 'Cell', 'Fixed', 'Pose old', 'Moments', 'No data']
    axes[2].set_yticks(range(len(METHODS)), short_labels)
    axes[2].set_xticks(range(5), ['Cell / 0', 'Cell / rand.', 'Cell / coh.', 'L2 / 0', 'L2 / coh.'], rotation=50, ha='right')
    axes[2].set_title('Minimum in-class analytic coverage')
    for i in range(len(METHODS)):
        for j in range(5):
            axes[2].text(j, i, f'{matrix[i,j]:.2f}', ha='center', va='center', fontsize=7,
                         color='black' if matrix[i,j] > .65 else 'white')
    fig.suptitle('All frozen designs: blue 10028, orange 10049, green 10076; circle average, triangle contrast', fontsize=8)
    fig.tight_layout()
    for extension in ['png', 'pdf']:
        fig.savefig(OUT/f'confirmation.{extension}', dpi=200, bbox_inches='tight',
                    **({'metadata': {'CreationDate': None, 'ModDate': None}} if extension == 'pdf' else {}))
    plt.close(fig)


def table(widths, setting_count, cfg):
    ncols = len(cfg['widths'])*len(cfg['angles'])
    caption = ('Frozen continuous-density and pose validation: relative half-width ranges over both targets and all four geometry replicates per stack. '
               f'All {setting_count} settings and all in-class signals are retained. The density, noise and pose limits are supplied assumptions.')
    scale_header = ' & '+' & '.join(r'\multicolumn{2}{c}{Scale '+f'{width:.2f}'+'}' for width in cfg['widths'])+r' \\'
    angle_header = 'EMPIAR & '+' & '.join(f'${angle:g}^\\circ$' for width in cfg['widths'] for angle in cfg['angles'])+r' \\'
    label = 'frozenuncertainty' if ARGS.study == 'continuous-v1' else 'frozenmoments'
    lines = [r'\begin{table}[t]', r'\centering\small',
             r'\caption{'+caption+'}', r'\label{tab:'+label+'}',
             r'\begin{tabular}{l'+'c'*ncols+'}', r'\toprule', scale_header, angle_header, r'\midrule']
    for dataset in ['10028', '10049', '10076']:
        cells = []
        for width in cfg['widths']:
            for angle in cfg['angles']:
                rows = [r for r in widths if r['dataset'] == dataset and r['width'] == width
                        and r['angle_degrees'] == angle and r['method'] == 'continuous_pose']
                lo = min(r['relative_half_width']['minimum'] for r in rows)
                hi = max(r['relative_half_width']['maximum'] for r in rows)
                cells.append(f'{lo:.2f}--{hi:.2f}')
        lines.append(dataset+' & '+' & '.join(cells)+r' \\')
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    filename = 'frozen-uncertainty-results.tex' if ARGS.study == 'continuous-v1' else 'frozen-moments-results.tex'
    (ROOT/'paper'/filename).write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
