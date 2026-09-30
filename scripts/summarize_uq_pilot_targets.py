#!/usr/bin/env python3
"""Summarize the entire locked feature family; refuse to omit missing cases."""
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
    lock = read_target_lock(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
    out = BASE/'pilot-selected-summary-v1'
    if (out/'summary.json').exists():
        raise RuntimeError('Preserve the completed summary')
    sources = {}; conditional = []; experimental = []; failures = []
    def read(path):
        raw = path.read_bytes(); row = json.loads(raw)
        if not row.get('complete') or row.get('error'):
            raise RuntimeError(f'Missing successful completion: {path}')
        sources[str(path.relative_to(ROOT))] = hashlib.sha256(raw).hexdigest()
        if row.get('reference_coverage_failure'):
            failures.append(str(path.relative_to(ROOT)))
        return row
    for dataset in lock['datasets']:
        ds = dataset['dataset']
        for feature in dataset['features']:
            name = feature['name']
            fp = BASE/'pilot-selected-fixed-v1'/f'{ds}-{name}.json'; fixed = read(fp)
            if fixed['target_lock_sha256'] != LOCK_SHA256:
                raise ValueError('Source target lock changed')
            target = next(r for r in fixed['targets'] if r['target'] == name)
            fit = target['fit']; no_data = 2*fit['target_norm']
            common = {'dataset': ds, 'feature': name, 'width_A': 20., 'band_endpoint_A': fixed['band_endpoint_A']}
            conditional.append(dict(common, pose_class='fixed_pose',
                relative_width=target['selected_relative_half_width'],
                minimum_reference_power=target['reference_check']['correct_sign_probability'],
                minimum_reference_coverage=target['reference_check']['analytic_coverage'],
                uses_no_data=target['reference_check']['uses_no_data']))
            for angle in [0, 1, 2]:
                ap = BASE/'pilot-selected-pose-v1'/f'{ds}-{name}-{angle}.json'; audit = read(ap)
                ep = BASE/'pilot-selected-enclosing-v1'/ap.name; enclosure = read(ep)
                if audit['source_fit_sha256'] != sources[str(fp.relative_to(ROOT))]:
                    raise ValueError('Source fit changed since pose audit')
                if enclosure['source_audit_sha256'] != sources[str(ap.relative_to(ROOT))]:
                    raise ValueError('Source audit changed since refinement')
                conditional.append(dict(common, pose_class='shift_only' if angle == 0 else f'pose_{angle}deg',
                    relative_width=enclosure['selected_relative_half_width'],
                    minimum_reference_power=enclosure['minimum_reference_power'],
                    minimum_reference_coverage=min(r['analytic_coverage'] for r in enclosure['reference_checks']),
                    uses_no_data=enclosure['uses_no_data']))
            applied = read(BASE/'pilot-selected-experimental-v1'/f'{ds}-{name}.json')
            if applied['fixed_fit_sha256'] != sources[str(fp.relative_to(ROOT))] or applied['target_lock_sha256'] != LOCK_SHA256:
                raise ValueError('Experimental source changed')
            if [r['pose_class'] for r in applied['records']] != ['fixed_pose','shift_only','pose_1deg','pose_2deg']:
                raise ValueError('All four experimental sensitivity classes are required')
            for r in applied['records']:
                experimental.append(dict(common, pose_class=r['pose_class'],
                    relative_width=r['relative_half_width'], center=r['interval_center'], half_width=r['interval_half_width'],
                    no_data_half_width=no_data, approximate_reference=r['approximate_reference_target'],
                    approximate_reference_inside=r['approximate_reference_inside_interval'],
                    excludes_zero=r['excludes_zero'], uses_no_data=r['uses_no_data'],
                    reference_inside_declared_class=applied['reference_inside_density_ball_numerical'],
                    inference_groups=applied['group_noise_calibration']['inference_groups'],
                    maximum_group_size=applied['group_noise_calibration']['maximum_group_size'],
                    calibrated_noise_sd=r['noise_sd_upper']))
    if len(conditional) != 48 or len(experimental) != 48:
        raise AssertionError('Incomplete feature family')
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in [('conditional', conditional), ('experimental', experimental)]:
        with (out/f'{name}.csv').open('w') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    labels = ['fixed_pose', 'shift_only', 'pose_1deg', 'pose_2deg']
    colors = ['#0072B2', '#009E73', '#E69F00', '#CC79A7']
    readable = ['Fixed poses', '0.5 Å shifts', '1° + 0.5 Å', '2° + 0.5 Å']
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), sharey=True)
    for ax, dataset in zip(axes, lock['datasets']):
        ds = dataset['dataset']
        for fi, feature in enumerate(dataset['features']):
            name = feature['name']; rows = [r for r in experimental if r['dataset'] == ds and r['feature'] == name]
            for pi, (label, color, text) in enumerate(zip(labels, colors, readable)):
                row = next(r for r in rows if r['pose_class'] == label); y = fi+(pi-1.5)*.14
                ax.errorbar(row['center']/row['no_data_half_width'], y,
                            xerr=row['half_width']/row['no_data_half_width'], fmt='o', color=color,
                            markersize=3, capsize=2, linewidth=1.1, label=text if fi == 0 else None)
            row = rows[0]; ref = row['approximate_reference']/row['no_data_half_width']
            ax.plot([ref, ref], [fi-.32, fi+.32], color='black', linestyle=':', linewidth=1.1,
                    label='Approx. map value' if fi == 0 else None)
        ax.set_title(f'EMPIAR-{ds}'); ax.axvline(0, color='.5', linewidth=.6); ax.grid(axis='x', alpha=.2)
        ax.set_xlabel('Feature / no-data half-width')
    axes[0].set_yticks(range(4), ['Pilot region 1', 'Pilot region 2', 'Pilot region 3', 'Matched center'])
    axes[0].invert_yaxis()
    handles, texts = axes[0].get_legend_handles_labels()
    fig.legend(handles, texts, loc='lower center', ncol=5, fontsize=8)
    fig.suptitle('Experimental conditional intervals; nuisance assumptions remain unverified', fontsize=11)
    fig.tight_layout(rect=(0, .09, 1, .94))
    for extension in ['pdf', 'png']:
        fig.savefig(out/f'experimental-intervals.{extension}', dpi=180)
    plt.close(fig)
    summary = {'complete': True, 'scope': 'Exploratory locked targets; known-noise calculations and conditional experimental intervals are different evidence',
               'target_lock_sha256': LOCK_SHA256, 'source_hashes': sources,
               'conditional': conditional, 'experimental': experimental,
               'reference_coverage_failure_records': failures,
               'conditional_family_alpha': .05, 'family_size_per_pose_class': 12,
               'experimental_budget': {'alpha_noise_each': (.045-1e-6)/12, 'beta_calibration_each': .005/12, 'delta_spectral_each': 1e-6/12},
               'interpretation': 'A processed map is not a coverage label. Targets have Gaussian sigma 20 A, not claimed reconstruction resolution. This is not a new frozen confirmation.'}
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print('COMPLETE', 'conditional', len(conditional), 'experimental', len(experimental), 'calculation_failures', len(failures), flush=True)


if __name__ == '__main__':
    main()
