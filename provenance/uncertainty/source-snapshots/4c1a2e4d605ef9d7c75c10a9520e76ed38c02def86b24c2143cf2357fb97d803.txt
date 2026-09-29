#!/usr/bin/env python3
"""Execute the frozen prediction-v1 protocol; never tune on the new exposures."""
import argparse
import csv
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import torch
from evaluate_uq_reconstruction import neural_prediction, metrics
from run_experiment import observations, predict
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT/'research/uncertainty/confirmation/prediction-v1'
SNAPSHOT = source_snapshot(ROOT, Path(__file__), ['scripts/evaluate_uq_reconstruction.py', 'scripts/run_experiment.py'])


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as file:
        for block in iter(lambda: file.read(16*1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def run(dataset):
    begin = time.perf_counter(); torch.set_num_threads(4)
    locked = json.loads((PROTOCOL/'locked-models.json').read_text())
    selection = json.loads((PROTOCOL/f'{dataset}-manifest.json').read_text())
    for path, digest in locked['models'][dataset]['files'].items():
        if sha(ROOT/path) != digest:
            raise AssertionError(f'Locked model drift: {path}')
    for path in ['scripts/evaluate_uq_reconstruction.py', 'scripts/run_experiment.py']:
        if sha(ROOT/path) != locked['source_hashes'][path]:
            raise AssertionError(f'Locked evaluation helper drift: {path}')
    csv_path = PROTOCOL/f'{dataset}-selection.csv'
    if sha(csv_path) != selection['selection_csv_sha256']:
        raise AssertionError('Frozen cohort selection drift')
    rows = list(csv.DictReader(csv_path.open())); groups = np.array([r['source_group'] for r in rows])
    previous = list(csv.DictReader((ROOT/'research/uncertainty/splits'/f'{dataset}.csv').open()))
    if set(groups) & {r['source_group'] for r in previous}:
        raise AssertionError('Exposure-group overlap')
    data = ROOT/'data/uncertainty/confirmation/prediction-v1'/dataset
    download = json.loads((data/'manifest.json').read_text())
    if download['selection_file_sha256'] != selection['indices_npy_sha256'] or download['count'] != 4096:
        raise AssertionError('Download does not match frozen selection')
    np.testing.assert_array_equal(np.load(data/'indices.npy'), [int(r['metadata_source_index']) for r in rows])
    metadata = np.load(data/'metadata.npz')
    np.testing.assert_array_equal(metadata['indices'], np.load(data/'indices.npy'))
    k, observed, transfer, _, _ = observations(data, 64, 2048, 42, True)
    if observed.shape != (4096, 1410):
        raise AssertionError('Unexpected observation count')
    classical = ROOT/'results/uncertainty/development/group-reconstruction'/dataset
    neural = ROOT/'results/uncertainty/development/neural-reconstruction'/dataset
    scale = locked['models'][dataset]['classical_global_scale']
    epoch = locked['epochs'][dataset]
    values = {}; errors = []; checks = []
    for method in ['gaussian', 'voxel', 'neural']:
        if method == 'neural':
            predicted = np.zeros_like(observed)
            for half in [0, 1]:
                p, _, check = neural_prediction(neural/f'half{half}', epoch, k, transfer)
                predicted += p/2; checks.append(check)
        else:
            coefficients = [np.load(classical/f'{method}-half{h}.npz') for h in [0, 1]]
            real = sum(x['real'] for x in coefficients)/2; imag = sum(x['imag'] for x in coefficients)/2
            predicted = predict(k.reshape(-1, 3), transfer.ravel(), real, imag,
                                argparse.Namespace(box=64, sigma=.5, radius=2), method).reshape(observed.shape)*scale
        values[method], error, power = metrics(predicted, observed); errors.append(error)
        print(dataset, method, values[method], flush=True)
    labels, inverse = np.unique(groups, return_inverse=True)
    grouped_error = np.stack([np.bincount(inverse, weights=e, minlength=len(labels)) for e in errors])
    grouped_power = np.bincount(inverse, weights=power, minlength=len(labels))
    counts = np.bincount(inverse, minlength=len(labels))
    seed = 609422+int(dataset); rng = np.random.default_rng(seed)
    choices = rng.integers(len(labels), size=(10000, len(labels)))
    ratios = grouped_error[:, choices].sum(axis=-1)/grouped_power[choices].sum(axis=-1)
    comparisons = {}
    for index, method in enumerate(['gaussian', 'voxel']):
        difference = ratios[index]-ratios[2]
        comparisons[method+'_minus_neural'] = {
            'estimate': values[method]['nmse']-values['neural']['nmse'],
            'percentile_95_marginal': np.quantile(difference, [.025, .975]).tolist(),
            'percentile_bonferroni_six_contrasts': np.quantile(difference, [.05/12, 1-.05/12]).tolist()}
    out = ROOT/'results/uncertainty/confirmation/prediction-v1'/dataset; out.mkdir(parents=True, exist_ok=True)
    with (out/'exposure-aggregates.csv').open('w') as file:
        writer = csv.writer(file, lineterminator='\n')
        writer.writerow(['source_group', 'particles', 'observed_power_sum', 'gaussian_error_sum', 'voxel_error_sum', 'neural_error_sum'])
        writer.writerows(zip(labels, counts, grouped_power, *grouped_error))
    np.savez(out/'particle-errors.npz', indices=np.load(data/'indices.npy'), errors=np.asarray(errors), power=power)
    result = {'stage': 'frozen additional-exposure prediction confirmation; conditional published poses/CTFs; not density coverage',
              'protocol_commit': '72a6dc0', 'protocol_sha256': sha(PROTOCOL/'PROTOCOL.md'),
              'locked_models_sha256': sha(PROTOCOL/'locked-models.json'), 'source_snapshot': SNAPSHOT,
              'dataset': dataset, 'particles': 4096, 'exposure_groups': len(labels), 'development_group_overlap': 0,
              'neural_epoch': epoch, 'metrics': values, 'paired_comparisons': comparisons, 'checkpoint_checks': checks,
              'bootstrap': {'resamples': 10000, 'seed': seed, 'unit': 'source exposure group', 'paired': True,
                            'bonferroni_contrasts': 6, 'coverage_claim': 'ordinary bootstrap approximation, not a finite-sample theorem'},
              'data_hashes': {name: sha(data/name) for name in ['images.npy', 'metadata.npz', 'indices.npy', 'manifest.json']},
              'selection_csv_sha256': sha(csv_path), 'seconds': time.perf_counter()-begin}
    (out/'metrics.json').write_text(json.dumps(result, indent=2)+'\n')
    # Public range provenance contains only identifiers and checksums, no pixels.
    provenance = ROOT/'provenance/uncertainty/confirmation/prediction-v1'; provenance.mkdir(parents=True, exist_ok=True)
    (provenance/f'{dataset}-download.json').write_text(json.dumps(download, indent=2)+'\n')
    print(dataset, 'FROZEN PREDICTION COMPLETE', result['seconds'], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--datasets', default='10028,10049,10076'); args = p.parse_args()
    for dataset in args.datasets.split(','):
        run(dataset)
