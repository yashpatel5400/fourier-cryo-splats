#!/usr/bin/env python3
"""Descriptive background diagnostics, not a learned noise certificate.

Uses the entire already-evaluated additional-exposure cohort. Background pixels
can include molecular signal, ice and preprocessing correlations. These summaries
cannot identify pure measurement noise or calibrate density confidence intervals.
"""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/uncertainty/development/background-diagnostics'
LAGS = [1, 2, 3, 4]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as file:
        for block in iter(lambda: file.read(16*1024**2), b''):
            h.update(block)
    return h.hexdigest()


def statistics(images):
    n, box, _ = images.shape
    axis = np.arange(box)-box//2
    y, x = np.meshgrid(axis, axis, indexing='ij')
    mask = x*x+y*y > (.43*box)**2
    centered = images-images[:, mask].mean(axis=1)[:, None, None]
    variance = np.mean(centered[:, mask]**2, axis=1)
    if np.any(variance <= 0):
        raise ValueError('Constant background image')
    normalized = centered/np.sqrt(variance)[:, None, None]
    values = [np.mean(normalized[:, mask]**4, axis=1)-3]
    columns = ['background_excess_kurtosis']
    for lag in LAGS:
        for direction, a, b, valid in [
            ('x', normalized[:, :, :-lag], normalized[:, :, lag:], mask[:, :-lag]&mask[:, lag:]),
            ('y', normalized[:, :-lag, :], normalized[:, lag:, :], mask[:-lag, :]&mask[lag:, :])]:
            values.append(np.mean((a*b)[:, valid], axis=1))
            columns.append(f'normalized_pair_product_{direction}_lag{lag}')
    return columns, np.stack(values, axis=1), np.sqrt(variance), int(mask.sum())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    snapshot = source_snapshot(ROOT, Path(__file__))
    summary = {'stage': 'post-hoc descriptive background diagnostics, not noise calibration',
               'source_snapshot': snapshot, 'datasets': [],
               'interpretation': ['The annular background can contain signal, ice and correlated preprocessing artifacts.',
                                  'Normalized pair products use each whole background variance, not a fitted stationary covariance.',
                                  'White-Gaussian controls pass through the same centering/normalization statistic, not the archive preprocessing.',
                                  'Exposure bootstrap intervals are approximate marginal intervals; no simultaneous or density-coverage claim.',
                                  'These diagnostics do not distinguish noise misspecification from residual signal.']}
    for dataset in ['10028', '10049', '10076']:
        base = ROOT/'data/uncertainty/confirmation/prediction-v1'/dataset
        images = np.load(base/'images.npy', mmap_mode='r')
        rows = list(csv.DictReader((ROOT/'research/uncertainty/confirmation/prediction-v1'/f'{dataset}-selection.csv').open()))
        groups = np.array([r['source_group'] for r in rows])
        if images.shape != (4096, 64, 64) or len(groups) != len(images):
            raise AssertionError('Unexpected diagnostic cohort')
        observed = []; controls = []; scales = []
        rng = np.random.default_rng(609551+int(dataset))
        for start in range(0, len(images), 128):
            block = np.asarray(images[start:start+128], dtype=float)
            columns, stats, scale, count = statistics(block)
            observed.append(stats); scales.append(scale)
            controls.append(statistics(rng.normal(size=block.shape))[1])
        observed = np.concatenate(observed); controls = np.concatenate(controls); scales = np.concatenate(scales)
        labels, inverse = np.unique(groups, return_inverse=True)
        counts = np.bincount(inverse)
        totals = np.stack([np.bincount(inverse, weights=observed[:, j]) for j in range(len(columns))])
        control_totals = np.stack([np.bincount(inverse, weights=controls[:, j]) for j in range(len(columns))])
        choices = rng.integers(len(labels), size=(10000, len(labels)))
        boot = totals[:, choices].sum(axis=-1)/counts[choices].sum(axis=-1)
        boot_control = control_totals[:, choices].sum(axis=-1)/counts[choices].sum(axis=-1)
        report = {'dataset': dataset, 'particles': len(images), 'exposure_groups': len(labels),
                  'background_pixels': count, 'seed': 609551+int(dataset), 'bootstrap_replicates': 10000,
                  'source_image_sha256': sha(base/'images.npy'),
                  'background_sd_percentiles': np.quantile(scales, [.01, .1, .5, .9, .99]).tolist(),
                  'statistics': []}
        for j, name in enumerate(columns):
            report['statistics'].append({'name': name, 'mean': float(observed[:, j].mean()),
                                        'median': float(np.median(observed[:, j])),
                                        'percentile95_exposure_bootstrap': np.quantile(boot[j], [.025, .975]).tolist(),
                                        'white_control_mean': float(controls[:, j].mean()),
                                        'white_control_percentile95': np.quantile(boot_control[j], [.025, .975]).tolist()})
        with (OUT/f'{dataset}-exposure-aggregates.csv').open('w') as file:
            writer = csv.writer(file, lineterminator='\n')
            writer.writerow(['source_group', 'particles']+[f'sum_{c}' for c in columns]+[f'control_sum_{c}' for c in columns])
            writer.writerows(zip(labels, counts, *totals, *control_totals))
        summary['datasets'].append(report)
        print(dataset, 'lag1 x/y', report['statistics'][1]['mean'], report['statistics'][2]['mean'], flush=True)
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    fig, axes = plt.subplots(1, 3, figsize=(9., 2.7), sharey=True)
    for ax, report in zip(axes, summary['datasets']):
        for direction, marker in [('x', 'o'), ('y', '^')]:
            records = [next(r for r in report['statistics'] if r['name'] == f'normalized_pair_product_{direction}_lag{lag}') for lag in LAGS]
            mean = np.array([r['mean'] for r in records]); interval = np.array([r['percentile95_exposure_bootstrap'] for r in records])
            ax.errorbar(LAGS, mean, yerr=np.stack([mean-interval[:, 0], interval[:, 1]-mean]),
                        marker=marker, label=direction, capsize=2, linewidth=1)
        ax.axhline(report['statistics'][1]['white_control_mean'], color='.5', linestyle=':', label='White control')
        ax.set_title('EMPIAR-'+report['dataset']); ax.set_xticks(LAGS); ax.set_xlabel('Lag (working pixels)')
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('Normalized background pair product'); axes[-1].legend(fontsize=7)
    fig.tight_layout()
    for extension in ['png', 'pdf']:
        fig.savefig(OUT/f'background.{extension}', dpi=200, bbox_inches='tight',
                    **({'metadata': {'CreationDate': None, 'ModDate': None}} if extension == 'pdf' else {}))
    plt.close(fig)


if __name__ == '__main__':
    main()
