#!/usr/bin/env python3
"""Compare every completed nonlinear stress candidate with both upper bounds."""
import json
import hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fourier_splats.uncertainty import bias_aware_half_width
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
OUT = BASE/'continuous-pose-adversaries'
rows = []; all_candidates = []; sources = {}; seconds = 0.
for dataset in ['10028', '10049', '10076']:
    path = OUT/f'{dataset}.json'; upper_path = BASE/'continuous-pose-moments'/f'{dataset}.json'
    data = json.loads(path.read_text()); upper = json.loads(upper_path.read_text())
    if len(data['records']) != 6 or any(len(r['candidates']) != 6 for r in data['records']):
        raise AssertionError('Every dataset, target, angle, sign and start must complete')
    for p in [path, upper_path]: sources[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    for record in data['records']:
        audit = next(r for r in upper['records'] if r['target'] == record['target'] and r['angle_degrees'] == record['angle_degrees'])
        seconds += record['seconds']; all_candidates += record['candidates']
        best = record['largest_feasible_bias_numerical']
        if best > audit['new_total_bias']+1e-6:
            raise AssertionError('Feasible nonlinear bias exceeds integrated-moment certificate')
        if any(c['maximum_pose_norm'] > 1+1e-12 for c in record['candidates']):
            raise AssertionError('Infeasible candidate pose')
        rows.append({'dataset':dataset, 'target':record['target'], 'angle_degrees':record['angle_degrees'],
                     'largest_feasible_bias_numerical':best,
                     'feasible_bias_over_spatial_maximum_upper':best/audit['old_total_bias'],
                     'feasible_bias_over_integrated_moment_upper':best/audit['new_total_bias'],
                     'improvement_over_coherent_pose':best/audit['coherent_pose_feasible_bias_numerical'],
                     'feasible_width_over_no_data':float(bias_aware_half_width(audit['noise_sd'],best)/audit['no_data_half_width']),
                     'integrated_upper_width_over_no_data':audit['new_uncapped_relative_half_width'],
                     'seconds':record['seconds']})
errors = [c['independent_field_relative_error'] for c in all_candidates if c['independent_field_relative_error'] is not None]
ratios = [r['feasible_bias_over_integrated_moment_upper'] for r in rows]
gains = [r['improvement_over_coherent_pose'] for r in rows]
summary = {'stage':'development feasible lower bounds, not global nonlinear optima',
           'sources':sources, 'analysis_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'settings':len(rows), 'optimized_candidates':len(all_candidates),
           'maximum_independent_field_relative_error':max(errors),
           'feasible_bias_over_integrated_upper_range':[min(ratios),max(ratios)],
           'gain_over_coherent_pose_range':[min(gains),max(gains)], 'total_seconds':seconds,
           'all_candidates_below_integrated_upper':True,
           'roundoff_note':'Float32 surrogate search; final pose projection and exact integral formulas use ordinary float64, not validated arithmetic.',
           'records':rows}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.9), sharex=True)
for dataset, color in zip(['10028', '10049', '10076'], ['#2166ac', '#d95f02', '#1b9e77']):
    for target, marker in [('center', 'o'), ('contrast', '^')]:
        selected = [r for r in rows if r['dataset'] == dataset and r['target'] == target]
        angles = [r['angle_degrees'] for r in selected]
        axes[0].plot(angles, [r['feasible_bias_over_integrated_moment_upper'] for r in selected],
                     marker=marker, color=color, linewidth=.9, markersize=4)
        axes[1].plot(angles, [r['improvement_over_coherent_pose'] for r in selected],
                     marker=marker, color=color, linewidth=.9, markersize=4)
axes[0].axhline(1, color='.5', linestyle=':', linewidth=.8)
axes[0].set_ylim(0, 1.06); axes[0].set_ylabel('Feasible bias / integrated upper bound')
axes[1].set_ylabel('Feasible bias / coherent-pose bias')
for ax in axes:
    ax.set_xticks([.5, 1, 2]); ax.set_xlabel('Rotation radius (degrees)')
fig.suptitle('Blue 10028, orange 10049, green 10076; circle average, triangle contrast', fontsize=8)
fig.tight_layout()
for extension in ['png', 'pdf']:
    fig.savefig(OUT/f'feasible-bias.{extension}', dpi=200, bbox_inches='tight',
                **({'metadata': {'CreationDate': None, 'ModDate': None}} if extension == 'pdf' else {}))
plt.close(fig)
print(json.dumps({k:v for k,v in summary.items() if k not in ['sources','records']},indent=2))
