#!/usr/bin/env python3
"""Re-audit all saved broad-feature development weights with integrated moments."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous_moments import integrated_cubic_remainder
from fourier_splats.uncertainty import bias_aware_half_width
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
OUT = BASE/'continuous-pose-moments'; OUT.mkdir(exist_ok=True)
snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py'])
rows = []
for dataset in ['10028', '10049', '10076']:
    fit_path = BASE/'continuous-quadrature-optimized'/f'{dataset}.json'
    fits = json.loads(fit_path.read_text()); cfg = fits['config']
    g = particle_geometry(ROOT, dataset, 'inference_half0', radius=cfg['frequency_radius'],
                          count=cfg['particles'], seed=cfg['seed'])
    op, pilot, _, noise = model(g, np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz'), 24)
    P = float(np.linalg.norm(pilot)); records = []; sources = {str(fit_path.relative_to(ROOT)): hashlib.sha256(fit_path.read_bytes()).hexdigest()}
    for name in ['continuous-pose', 'continuous-pose-large']:
        path = BASE/name/f'{dataset}.json'; data = json.loads(path.read_text())
        if len(data['records']) != (4 if name == 'continuous-pose' else 6):
            raise AssertionError('All original development cases must be complete')
        sources[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        for original in data['records']:
            target = original['target']; degrees = original['angle_degrees']
            fit = next(t for t in fits['targets'] if t['target'] == target and t['width_fraction_field'] == .07)
            saved = np.load(BASE/'continuous-quadrature-optimized'/f'{dataset}-{target}-0.07-weights.npz')
            np.testing.assert_array_equal(saved['indices'], g['indices'])
            w = saved['weights']; no_data = 2*fit['fit']['target_norm']
            refined = integrated_cubic_remainder(g['k'], g['q'], g['ctf'], w, noise, np.deg2rad(degrees), .01/24, 2., P)
            refined['per_particle_operator_third_bound'] = refined['per_particle_operator_third_bound'].tolist()
            bias = original['audit']['density_bias']+original['audit']['pose_polynomial_bias']+refined['remainder_bias']
            sd = float(np.linalg.norm(w)); half = float(bias_aware_half_width(sd, bias))
            lower = original['checks'][-1]['continuous_density_ball_eliminated_bias']
            if lower > bias+1e-6:
                raise AssertionError('Feasible continuous bias exceeds refined upper bound')
            if refined['remainder_bias'] > original['audit']['remainder_bias']+1e-8:
                raise AssertionError('Moment envelope should improve the spatial-maximum envelope')
            record = {'dataset': dataset, 'target': target, 'angle_degrees': degrees,
                      'old_total_bias': original['audit']['total_bias'], 'new_total_bias': bias,
                      'old_remainder_bias': original['audit']['remainder_bias'], 'integrated_moment_bound': refined,
                      'no_data_half_width': no_data, 'noise_sd': sd,
                      'old_uncapped_relative_half_width': original['half_width_before_fallback']/no_data,
                      'new_uncapped_relative_half_width': half/no_data,
                      'new_selected_relative_half_width': min(half/no_data, 1.),
                      'new_uses_no_data': half >= no_data,
                      'coherent_pose_feasible_bias_numerical': lower,
                      'feasible_bias_over_new_upper': lower/bias,
                      'reference_coverages_before_fallback': [float(norm.cdf((half-c['reference_bias'])/sd)
                                                                       -norm.cdf((-half-c['reference_bias'])/sd)) for c in original['checks']],
                      'weights_sha256': hashlib.sha256(w.tobytes()).hexdigest()}
            records.append(record); rows.append(record)
    (OUT/f'{dataset}.json').write_text(json.dumps({'stage': 'post-hoc development refinement; frozen continuous-v1 unchanged',
        'source_snapshot': snapshot, 'sources': sources, 'records': records}, indent=2)+'\n')

plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True)
for target, ax in zip(['center', 'contrast'], axes):
    for dataset, color in zip(['10028', '10049', '10076'], ['#2166ac', '#d95f02', '#1b9e77']):
        part = sorted([r for r in rows if r['target'] == target and r['dataset'] == dataset], key=lambda r:r['angle_degrees'])
        x = [r['angle_degrees'] for r in part]
        ax.plot(x, [r['new_uncapped_relative_half_width'] for r in part], '-o', color=color, markersize=3, label=dataset)
        ax.plot(x, [r['old_uncapped_relative_half_width'] for r in part], ':', color=color, alpha=.7)
    ax.axhline(1, color='.35', linestyle='--', linewidth=.8)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xticks([.1, .5, 1, 2, 5], ['0.1', '0.5', '1', '2', '5'])
    ax.set_xlabel('Rotation radius (degrees)'); ax.set_title('Central average' if target == 'center' else 'Axial contrast')
axes[0].set_ylabel('Uncapped half-width / no-data half-width')
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=3, frameon=False,
           title='Solid: integrated moments; dotted: spatial maximum; dashed: no data')
fig.tight_layout(rect=[0, 0, 1, .81])
fig.savefig(OUT/'moment-remainder.png', dpi=200, bbox_inches='tight')
fig.savefig(OUT/'moment-remainder.pdf', bbox_inches='tight', metadata={'CreationDate':None, 'ModDate':None})
summary = {'settings': len(rows), 'old_no_data_count': sum(r['old_uncapped_relative_half_width'] >= 1 for r in rows),
           'new_no_data_count': sum(r['new_uses_no_data'] for r in rows),
           'remainder_reduction_range': [min(r['old_remainder_bias']/r['integrated_moment_bound']['remainder_bias'] for r in rows),
                                         max(r['old_remainder_bias']/r['integrated_moment_bound']['remainder_bias'] for r in rows)],
           'one_degree_relative_width_range': [min(r['new_selected_relative_half_width'] for r in rows if r['angle_degrees']==1),
                                               max(r['new_selected_relative_half_width'] for r in rows if r['angle_degrees']==1)]}
(OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
