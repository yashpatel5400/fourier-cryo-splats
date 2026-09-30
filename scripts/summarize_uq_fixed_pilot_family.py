#!/usr/bin/env python3
"""Summarize the entire completed fixed-pose family, separately from pose audits."""
import hashlib
import json
from pathlib import Path
from fourier_splats.uq_pilot_targets import read_target_lock, LOCK_SHA256

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'


def main():
    lock = read_target_lock(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
    out = BASE/'pilot-selected-fixed-summary-v1'; out.mkdir(parents=True, exist_ok=True)
    path = out/'summary.json'
    if path.exists():
        raise RuntimeError('Preserve existing summary')
    rows = []; sources = {}
    for dataset in lock['datasets']:
        for feature in dataset['features']:
            fp = BASE/'pilot-selected-fixed-v1'/f'{dataset["dataset"]}-{feature["name"]}.json'
            result = json.loads(fp.read_text())
            if not result.get('complete') or result.get('error') or result['target_lock_sha256'] != LOCK_SHA256:
                raise ValueError('Every locked fixed-pose outcome must be complete')
            t = next(r for r in result['targets'] if r['target'] == feature['name']); f = t['fit']; c = t['reference_check']
            rows.append({'dataset': dataset['dataset'], 'feature': feature['name'],
                'relative_half_width': t['selected_relative_half_width'],
                'reference_sign_power': c['correct_sign_probability'], 'reference_analytic_coverage': c['analytic_coverage'],
                'uses_no_data': c['uses_no_data'], 'optimizer_converged': f['converged'],
                'relative_sum_gap': f['history'][-1]['relative_gap'], 'seconds': result['seconds']})
            sources[str(fp.relative_to(ROOT))] = hashlib.sha256(fp.read_bytes()).hexdigest()
    if len(rows) != 12:
        raise ValueError('Incomplete target family')
    grouped = []
    for dataset in lock['datasets']:
        ds = dataset['dataset']; group = [r for r in rows if r['dataset'] == ds]
        grouped.append({'dataset': ds, 'features': len(group),
            'minimum_relative_half_width': min(r['relative_half_width'] for r in group),
            'maximum_relative_half_width': max(r['relative_half_width'] for r in group),
            'minimum_reference_sign_power': min(r['reference_sign_power'] for r in group),
            'maximum_relative_sum_gap': max(r['relative_sum_gap'] for r in group),
            'minimum_fit_seconds': min(r['seconds'] for r in group), 'maximum_fit_seconds': max(r['seconds'] for r in group),
            'converged': sum(r['optimizer_converged'] for r in group)})
    result = {'complete': True, 'target_lock_sha256': LOCK_SHA256, 'records': rows, 'by_dataset': grouped,
        'source_hashes': sources, 'scope': 'All twelve fixed-pose, known-noise conditional feature calculations; pose family and fresh noise applications remain separate.',
        'width_A_sigma': 20., 'particles': 128, 'frequency_radius': 12, 'alpha_total_each': .05/12}
    path.write_text(json.dumps(result, indent=2)+'\n')
    lines = [r'\begin{table}[t]', r'\centering', r'\setlength{\tabcolsep}{4.5pt}',
        r'\caption{Completed fixed-pose, known-noise calculations for all twelve locked targets (four per stack). Width is relative to no data; power uses the processed-map generator. These are conditional calculations, not experimental coverage.}',
        r'\label{tab:lockedfixed}', r'\begin{tabular}{lrrr}', r'\toprule',
        r'Stack & Width range & Min. sign power & Max. gap\\', r'\midrule']
    for r in grouped:
        power = '$>0.99999$' if r['minimum_reference_sign_power'] > .99999 else f'{r["minimum_reference_sign_power"]:.4f}'
        lines.append(f'{r["dataset"]} & {r["minimum_relative_half_width"]:.4f}--{r["maximum_relative_half_width"]:.4f} & {power} & {r["maximum_relative_sum_gap"]:.5f}'+r'\\')
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    (ROOT/'paper/tables/pilot-fixed-family.tex').write_text('\n'.join(lines)+'\n')
    print(json.dumps(grouped, indent=2))


if __name__ == '__main__':
    main()
