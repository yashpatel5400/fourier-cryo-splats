#!/usr/bin/env python3
"""Display every completed cubic design in simulation and on reused images.

No fit, estimator selection, new calibration, or new coverage experiment.
"""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
FAMILIES = ['cubic-weight-probe', 'cubic-coordinate-probe', 'cubic-subspace-probe',
            'cubic-enrichment-probe', 'joint-cubic-design-probe']
LABELS = ['Original cubic', 'Coordinate metric', 'Fixed subspace', 'Enriched subspace',
          'Joint density–pose']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-stem', type=Path, default=ROOT/'paper/figures/cubic-design-comparison')
    args = p.parse_args()
    stem = args.output_stem
    if any(stem.with_suffix(suffix).exists() for suffix in ['.json', '.pdf', '.png']):
        raise FileExistsError('Preserve previous figure and input manifest')
    inputs = {}

    def read(relative):
        path = BASE/relative
        r = json.loads(path.read_text())
        if not r.get('complete') or r.get('error') or r.get('scientific_run_complete') is False:
            raise ValueError(f'Incomplete or failed input: {relative}')
        inputs[str(path.relative_to(ROOT))] = sha(path)
        return r

    old = read('cubic-experimental-application-v1/10049.json')['records']
    centered = read('centered-noise-calibration-v1/10049.json')['cubic']
    newer = read('joint-enriched-experimental-v1/10049.json')['records']
    registered = read('registered-target-sensitivity-v1/10049.json')['cubic']
    rows = []
    for i, family in enumerate(FAMILIES):
        fit = read(f'{family}/10049-pilot_region_1-1.json')
        if i < 3:
            raw = next(r for r in old if r['estimator'] == family)
            cen = next(r for r in centered if r['estimator'] == family)
            power = min(r['correct_sign_probability'] for r in registered if r['estimator'] == family)
            old_power = fit['minimum_reference_power']
        else:
            raw = next(r for r in newer if r['estimator'] == family and r['calibration_method'] == 'uncentered')
            cen = next(r for r in newer if r['estimator'] == family and r['calibration_method'] == 'centered')
            power = fit['minimum_reference_power_by_frame']['registered']
            old_power = fit['minimum_reference_power_by_frame']['original']
        rows.append(dict(estimator=family, label=LABELS[i],
            simulation_minimum_power_original_frame=old_power,
            simulation_minimum_power_registered_frame=power,
            experimental_uncentered=dict(center=raw['interval_center'], half_width=raw['interval_half_width'],
                excludes_zero=raw['excludes_zero'], uses_no_data=raw['uses_no_data']),
            experimental_centered=dict(center=cen['interval_center'], half_width=cen['interval_half_width'],
                excludes_zero=cen['excludes_zero'], uses_no_data=cen['uses_no_data']),
            registered_reference_target=raw['registered_reference_target']))
    references = [r['registered_reference_target'] for r in rows]
    np.testing.assert_allclose(references, references[0], rtol=1e-12)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.3, 4.0), gridspec_kw={'width_ratios': [1, 1.45]})
    y = np.arange(len(rows))
    for key, marker, color, label in [
        ('simulation_minimum_power_original_frame', 'x', '#7c7c7c', 'Original reference frame'),
        ('simulation_minimum_power_registered_frame', 'o', '#147d92', 'Pilot-registered frame')]:
        ax.scatter([r[key] for r in rows], y, s=34, marker=marker, color=color, label=label, zorder=3)
    ax.set(yticks=y, yticklabels=LABELS, xlim=(-.04, 1.04), xlabel='Minimum simulated sign power',
           title='Supplied simulation noise')
    ax.invert_yaxis(); ax.grid(axis='x', alpha=.2)
    ax.legend(loc='upper center', bbox_to_anchor=(.5, -.24), frameon=False, fontsize=9)
    for key, offset, color, label in [
        ('experimental_uncentered', -.13, '#3066be', 'Raw calibration'),
        ('experimental_centered', .13, '#cf7023', 'Centered calibration')]:
        bx.errorbar([r[key]['center'] for r in rows], y+offset,
                    xerr=[r[key]['half_width'] for r in rows], fmt='o', ms=4,
                    color=color, capsize=3, lw=1.5, label=label)
    bx.axvline(0, color='black', ls=':', lw=1)
    bx.axvline(references[0], color='#147d92', ls='--', lw=1, label='Registered approximate reference')
    bx.set(yticks=y, yticklabels=[], xlabel='Feature value and conditional interval',
           title='Reused experimental images')
    bx.invert_yaxis(); bx.grid(axis='x', alpha=.2)
    bx.legend(loc='upper center', bbox_to_anchor=(.5, -.24), frameon=False, fontsize=9)
    fig.subplots_adjust(left=.17, right=.99, top=.90, bottom=.36, wspace=.16)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ['.pdf', '.png']:
        fig.savefig(stem.with_suffix(suffix), dpi=200)
    plt.close(fig)
    manifest = dict(complete=True, input_hashes=inputs, plotting_source_sha256=sha(Path(__file__)), rows=rows,
        scope='Descriptive plot of all five completed designs; supplied-noise power is not experimental calibration. The same reused observations and original per-procedure error allocations are retained; no minimum across calibration procedures is selected.',
        outputs={suffix:sha(stem.with_suffix(suffix)) for suffix in ['.pdf', '.png']})
    stem.with_suffix('.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(dict(figure=str(stem), rows=len(rows))))


if __name__ == '__main__':
    main()
