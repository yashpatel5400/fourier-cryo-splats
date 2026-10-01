#!/usr/bin/env python3
"""Replay saved cone witnesses without rerunning the optimizer."""
import hashlib
import json
from pathlib import Path
import numpy as np
from fourier_splats.uq_power_cone_witness import power_witness_log_upper

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results/uncertainty/development/paired-power-enlarged-cone-v1'


def main():
    summary_path = BASE / 'summary.json'
    summary = json.loads(summary_path.read_text())
    assert summary['complete'] and len(summary['cases']) == 36
    rows = []
    for dataset, array_record in summary['arrays'].items():
        path = ROOT / array_record['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == array_record['sha256']
        with np.load(path) as arrays:
            powers = np.abs(arrays['fourier']) ** 2
            signal = arrays['signal']
            np.testing.assert_allclose(signal, powers[0, :64].mean(0), rtol=1e-9, atol=1e-12)
            for case in summary['cases']:
                if case['dataset'] != dataset:
                    continue
                assert 'error' not in case
                name = case['candidate']; n = case['total_views']
                candidate = ['true_map', 'region_half_removed', 'region_removed'].index(name)
                key = f"{name}_{case['added_views']}"
                weights = arrays[key + '_coefficients']
                assert weights.shape == (n,) and weights.min() >= 0
                approximation = np.einsum('ij,i->j', powers[candidate, :n], weights)
                np.testing.assert_allclose(approximation, arrays[key + '_approximation'], rtol=1e-12, atol=1e-12)
                residual = np.max(np.abs(approximation-signal)/np.maximum(signal, 1e-6*signal.max()))
                np.testing.assert_allclose(residual, case['maximum_scaled_residual'], rtol=1e-8, atol=2e-13)
                discrepancies = []
                for variance in [1., 2., 4.]:
                    values = [power_witness_log_upper(signal, approximation, transfer, variance)
                              for transfer in arrays['transfer_squared']]
                    discrepancies.append(abs(16*sum(values)-case['expected_log_128_upper'][str(variance)]))
                rows.append(dict(dataset=dataset, candidate=name, views=n,
                    minimum_coefficient=float(weights.min()), maximum_scaled_residual=float(residual),
                    maximum_expected_log_replay_difference=max(discrepancies)))
    record = dict(scope='Independent saved-coefficient replay, not validated floating-point arithmetic.',
        source_summary_sha256=hashlib.sha256(summary_path.read_bytes()).hexdigest(),
        cases=rows, complete=True)
    out = ROOT / 'provenance/uncertainty/enlarged-power-cone-verification.json'
    if out.exists():
        raise RuntimeError('Preserve prior verification')
    out.write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(dict(cases=len(rows), maximum_replay_difference=max(
        r['maximum_expected_log_replay_difference'] for r in rows))))


if __name__ == '__main__':
    main()
