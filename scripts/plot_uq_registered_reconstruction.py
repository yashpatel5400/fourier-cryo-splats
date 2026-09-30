#!/usr/bin/env python3
"""Plot all three stacks, four methods, and both reference-frame FSCs."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from fourier_splats.fsc import resolution

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development/reference-registered-comparison-v1'
METHODS = [('gaussian', 'Gaussian; supplied poses', '#267d9c'),
           ('voxel', 'Voxel; supplied poses', '#b66d19'),
           ('neural', 'Neural; supplied poses', '#7a57a5'),
           ('relion', 'RELION; unknown poses', '#327a42')]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-stem', type=Path, default=ROOT/'paper/figures/registered-reconstruction-fsc')
    args = parser.parse_args()
    stem = args.output_stem
    if any(stem.with_suffix(s).exists() for s in ['.pdf', '.png', '.json']):
        raise FileExistsError('Preserve existing figure and source manifest')
    cases = []
    hashes = {}
    for ds in ['10028', '10049', '10076']:
        for method, _, _ in METHODS:
            path = BASE/f'{ds}-{method}.json'
            record = json.loads(path.read_text())
            if not record.get('complete') or record.get('error') or not record.get('scientific_run_complete'):
                raise ValueError(f'Incomplete or failed evaluation: {path}')
            hashes[str(path.relative_to(ROOT))] = sha(path)
            ap = ROOT/record['arrays_file']
            digest = sha(ap)
            if digest != record['arrays_sha256']:
                raise ValueError(f'Changed curve array: {ap}')
            hashes[str(ap.relative_to(ROOT))] = digest
            arrays = np.load(ap)
            curves = {frame: arrays[frame] for frame in ['original', 'registered']}
            for frame, curve in curves.items():
                if curve.shape != (30, 4) or not np.isfinite(curve).all():
                    raise ValueError('Expected all 30 finite original shells')
                metrics = record['fsc'][frame]
                np.testing.assert_allclose(np.mean(curve[:, 2]), metrics['mean_fsc'], rtol=1e-12, atol=1e-12)
                for threshold, key in [(.143, 'threshold_0143'), (.5, 'threshold_05')]:
                    checked = resolution(curve, threshold)
                    if checked['status'] != metrics[key]['status']:
                        raise ValueError('FSC threshold status changed')
                    if checked['angstrom'] is not None:
                        np.testing.assert_allclose(checked['angstrom'], metrics[key]['angstrom'], rtol=1e-12)
            cases.append(dict(dataset=ds, method=method, curves=curves,
                              metrics=record['fsc'], relion_converged=record.get('relion_converged')))
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    lower_limit = min(-.2, min(np.min(curve[:, 2]) for c in cases for curve in c['curves'].values())-.03)
    fig, axes = plt.subplots(2, 3, figsize=(8.8, 4.6), sharey=True)
    for col, ds in enumerate(['10028', '10049', '10076']):
        for row, frame in enumerate(['original', 'registered']):
            ax = axes[row, col]
            for method, label, color in METHODS:
                case = next(c for c in cases if c['dataset'] == ds and c['method'] == method)
                curve = case['curves'][frame]
                ax.plot(curve[:, 1], curve[:, 2], color=color, lw=1.5,
                        ls='--' if method == 'relion' else '-', label=label)
            ax.axhline(.143, color='#555555', lw=.7, ls=':')
            ax.axhline(.5, color='#aaaaaa', lw=.6, ls=':')
            ax.set(ylim=(lower_limit, 1.04), xlim=(0, curve[-1, 1]))
            ax.grid(alpha=.14)
            if row == 0:
                ax.set_title(f'EMPIAR-{ds}')
            else:
                ax.set_xlabel('Spatial frequency (Å⁻¹)')
            if col == 0:
                ax.set_ylabel(('Original reference' if row == 0 else 'Pilot-registered reference')+'\nFSC')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, bbox_to_anchor=(.53, .005))
    fig.subplots_adjust(left=.10, right=.99, top=.91, bottom=.21, hspace=.19, wspace=.12)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ['.pdf', '.png']:
        fig.savefig(stem.with_suffix(suffix), dpi=200)
    plt.close(fig)
    rows = [{k:v for k,v in case.items() if k != 'curves'} for case in cases]
    manifest = dict(complete=True, input_hashes=hashes, plotting_source_sha256=sha(Path(__file__)),
        cases=rows, method_stack_cases=len(cases), displayed_curves=2*len(cases),
        scope='All original and registered reference curves. This is agreement after declared pilot-only registration, not independent reconstruction accuracy or density coverage. Reference 10076 is a Class A assembly state; supplied-pose and unknown-pose tasks differ. Censored thresholds are not measured resolutions.',
        outputs={suffix:sha(stem.with_suffix(suffix)) for suffix in ['.pdf', '.png']})
    stem.with_suffix('.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(dict(figure=str(stem), method_stack_cases=len(cases), curves=2*len(cases))))


if __name__ == '__main__':
    main()
