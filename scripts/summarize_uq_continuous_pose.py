#!/usr/bin/env python3
"""Plot all completed continuous nonlinear-pose development audits."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fourier_splats.uncertainty import bias_aware_half_width
ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'
OUT = BASE/'continuous-summary'; OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
summary = []
for dataset, color in zip(['10028', '10049', '10076'], ['#2166ac', '#d95f02', '#1b9e77']):
    source = json.loads((BASE/'continuous-pose'/f'{dataset}.json').read_text())
    fixed = json.loads((BASE/'continuous-quadrature-optimized'/f'{dataset}.json').read_text())
    assert len(source['records']) == 4
    for name, ax in zip(['center', 'contrast'], axes):
        nominal = next(r for r in fixed['targets'] if r['target'] == name and r['width_fraction_field'] == .07)
        records = sorted([r for r in source['records'] if r['target'] == name], key=lambda r: r['angle_degrees'])
        no_data = 2*nominal['fit']['target_norm']; upper = [nominal['width_fraction_of_no_data']]; lower = [upper[0]]
        for row in records:
            check = row['checks'][-1]
            upper.append(row['width_fraction_of_no_data'])
            lower.append(float(bias_aware_half_width(row['noise_sd'], check['continuous_density_ball_eliminated_bias'])/no_data))
            summary.append({'dataset': dataset, 'target': name, 'angle': row['angle_degrees'], 'relative_width': upper[-1],
                            'feasible_lower_relative_width': lower[-1], 'seconds': row['seconds'],
                            'minimum_reference_coverage': min(c['reference_coverage_before_fallback'] for c in row['checks'])})
        ax.plot([0, .1, .5], upper, '-o', color=color, markersize=3, label=dataset)
        ax.plot([0, .1, .5], lower, '--o', color=color, markersize=2)
        ax.set_title('Central average' if name == 'center' else 'Axial contrast'); ax.set_xlabel('Rotation radius (degrees)')
        ax.set_ylim(0, .52); ax.set_xticks([0, .1, .5])
axes[0].set_ylabel('Half-width / no-data half-width')
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=3, frameon=False,
           title='Solid: uniform upper bound; dashed: one feasible-pose lower bound')
fig.tight_layout(rect=[0, 0, 1, .79])
fig.savefig(OUT/'continuous-pose.png', dpi=200, bbox_inches='tight')
fig.savefig(OUT/'continuous-pose.pdf', bbox_inches='tight', metadata={'CreationDate': None, 'ModDate': None})
(OUT/'continuous-pose-summary.json').write_text(json.dumps(summary, indent=2)+'\n')
