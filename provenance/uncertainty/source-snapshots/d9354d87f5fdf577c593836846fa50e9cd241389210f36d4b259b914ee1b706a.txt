#!/usr/bin/env python3
"""Recalibrate every locked interval on the reserved, previously unused exposures."""
import csv
import json
from pathlib import Path
import time
import numpy as np
from fourier_splats.physics import fft_center
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_group_noise import grouped_estimator_variance_upper
from fourier_splats.uq_noise_calibration import uncenter_fourier_weights
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_provenance import source_snapshot
from freeze_uq_noise_models import ROOT, COHORT, BASE, sha, verify_lock


def completed(path):
    row = json.loads(path.read_text())
    if not row.get('complete') or row.get('error'):
        raise ValueError(f'Completed outcome required: {path}')
    return row


def main():
    lock = verify_lock(); model_hash = sha(COHORT/'locked-models.json')
    # Check the complete planned family before creating any outcome record.
    for estimator in lock['estimators']:
        completed(BASE/f'pilot-selected-experimental-v1/{estimator["dataset"]}-{estimator["feature"]}.json')
    for dataset in ['10028', '10049', '10076']:
        completed(ROOT/f'provenance/uncertainty/noise-calibration-v1/{dataset}-download.json')
    out = ROOT/'results/uncertainty/confirmation/noise-calibration-v1'
    out.mkdir(parents=True, exist_ok=True)
    if any(out.glob('*.json')):
        raise RuntimeError('Preserve previous outcomes; this runner does not overwrite partial runs')
    snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/freeze_uq_noise_models.py',
        'research/uncertainty/confirmation/noise-calibration-v1/PROTOCOL.md'])
    all_rows = []
    for dataset in ['10028', '10049', '10076']:
        start = time.perf_counter()
        dp = ROOT/f'provenance/uncertainty/noise-calibration-v1/{dataset}-download.json'
        provenance = completed(dp)
        if provenance['model_lock_sha256'] != model_hash:
            raise ValueError('Calibration download does not follow this lock')
        for name, expected in provenance['cached_files'].items():
            if sha(ROOT/name) != expected:
                raise ValueError(f'Cached download changed: {name}')
        selection = list(csv.DictReader((COHORT/f'{dataset}-selection.csv').open()))
        labels = [r['source_group'] for r in selection]
        old = list(csv.DictReader((ROOT/f'research/uncertainty/splits/{dataset}.csv').open()))
        previous = list(csv.DictReader((ROOT/f'research/uncertainty/confirmation/prediction-v1/{dataset}-selection.csv').open()))
        excluded = {r['source_group'] for r in old+previous}
        if len(labels) != 128 or len(set(labels)) != 128 or set(labels) & excluded:
            raise ValueError('Fresh exposure disjointness failed')
        g = particle_geometry(ROOT, dataset, 'inference_half0', radius=12, count=128, seed=609315)
        original = completed(BASE/f'experimental-noise-grouped/{dataset}.json')
        amplitude = original['pilot_amplitude_raw_units']
        images = np.load(ROOT/f'data/uncertainty/confirmation/noise-calibration-v1/{dataset}/images.npy')
        q = g['q'][0].astype(int); center = images.shape[-1]//2
        ft = fft_center(images.astype(float))[:, q[:, 1]+center, q[:, 0]+center]
        calibration = np.concatenate([ft.real, ft.imag], axis=1)/amplitude
        phases = g['manifest']['data_sign']*np.exp(-2j*np.pi*np.einsum('nqi,ni->nq', g['q'], g['translations']))
        result = {'complete': False, 'dataset': dataset, 'model_lock_sha256': model_hash,
            'target_lock_sha256': lock['target_lock_sha256'], 'source_snapshot': snapshot,
            'download_provenance_sha256': sha(dp), 'calibration_groups': labels,
            'calibration_group_count': 128, 'excluded_group_overlap': 0,
            'scope': 'Fresh noise calibration of previously fixed exploratory estimators, not fresh density coverage validation',
            'family_size_per_pose_class': 12, 'alpha_noise': (.045-1e-6)/12,
            'beta_calibration': .005/12, 'delta_spectral': 1e-6/12,
            'unverified_conditions': original['unverified_conditions']+[
                'Arbitrary jointly Gaussian within-exposure dependence, independent exposures and one common marginal raw-coordinate covariance.',
                'A shared density remains assumed; the actual 10076 particles are heterogeneous.',
                'The inference images were previously inspected; only the calibration cohort is new.'],
            'features': [], 'seconds': None}
        path = out/f'{dataset}.json'
        def save():
            path.write_text(json.dumps(result, indent=2)+'\n')
        save()
        try:
            for estimator in [e for e in lock['estimators'] if e['dataset'] == dataset]:
                name = estimator['feature']; fp = ROOT/estimator['fit']; wp = ROOT/estimator['weights']
                fixed = completed(fp); fitted = next(r for r in fixed['targets'] if r['target'] == name)
                ap = BASE/f'pilot-selected-experimental-v1/{dataset}-{name}.json'; applied = completed(ap)
                if (applied['fixed_fit_sha256'] != sha(fp) or applied['source_weights_sha256'] != sha(wp)
                        or applied['target_lock_sha256'] != lock['target_lock_sha256']):
                    raise ValueError('Original application does not use the locked estimator')
                for script in ['scripts/apply_uq_pilot_targets.py', 'src/fourier_splats/uq_group_noise.py']:
                    if applied['source_snapshot']['sources'][script]['sha256'] != lock['files'][script]:
                        raise ValueError('Original application used a different calculation')
                saved = np.load(wp); np.testing.assert_array_equal(saved['indices'], g['indices'])
                raw = uncenter_fourier_weights(saved['weights'].reshape(128, -1), phases, float(saved['noise_std']))
                calibrated = grouped_estimator_variance_upper(calibration, raw, g['groups'], .005/12)
                sd = calibrated['noise_sd_upper']; no_data = 2*fitted['fit']['target_norm']
                pilot_target = fitted['pilot_expected_feature']
                row = {'target': name, 'source_application_sha256': sha(ap),
                    'variance_calibration': calibrated, 'old_noise_sd_upper': applied['group_noise_calibration']['noise_sd_upper'],
                    'fresh_to_old_sd_ratio': sd/applied['group_noise_calibration']['noise_sd_upper'],
                    'reference_pilot_L2_distance': applied['reference_pilot_L2_distance'],
                    'reference_inside_density_ball_numerical': applied['reference_inside_density_ball_numerical'], 'records': []}
                if [r['pose_class'] for r in applied['records']] != ['fixed_pose', 'shift_only', 'pose_1deg', 'pose_2deg']:
                    raise ValueError('All four sensitivity classes required')
                for record in applied['records']:
                    source = ROOT/record['bias_source']
                    if sha(source) != record['bias_source_sha256']:
                        raise ValueError('Source bias certificate changed')
                    bias = record['bias_upper']; half = bias_aware_half_width_stable(sd, bias, result['alpha_noise'])
                    fallback = bool(half >= no_data); half = min(half, no_data)
                    center = pilot_target if fallback else record['raw_observed_center']
                    fresh = dict(record, noise_sd_upper=sd, interval_center=center,
                        interval_half_width=half, relative_half_width=half/no_data,
                        raw_half_width=bias_aware_half_width_stable(sd, bias, result['alpha_noise']),
                        uses_no_data=fallback, excludes_zero=bool(abs(center) > half),
                        approximate_reference_inside_interval=bool(abs(record['approximate_reference_target']-center) <= half))
                    row['records'].append(fresh)
                    all_rows.append(dict(dataset=dataset, feature=name, pose_class=fresh['pose_class'],
                        fresh_to_old_sd_ratio=row['fresh_to_old_sd_ratio'], noise_sd_upper=sd,
                        relative_half_width=fresh['relative_half_width'], excludes_zero=fresh['excludes_zero'],
                        uses_no_data=fallback, approximate_reference_inside=fresh['approximate_reference_inside_interval']))
                result['features'].append(row); save()
                print('DONE', dataset, name, 'fresh/old SD', row['fresh_to_old_sd_ratio'], flush=True)
            result.update(complete=True, seconds=time.perf_counter()-start); save()
        except Exception as exc:
            result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise
    if len(all_rows) != 48:
        raise AssertionError('Incomplete family')
    with (out/'summary.csv').open('w') as file:
        writer = csv.DictWriter(file, fieldnames=all_rows[0], lineterminator='\n')
        writer.writeheader(); writer.writerows(all_rows)
    (out/'summary.json').write_text(json.dumps({'complete': True, 'records': all_rows,
        'model_lock_sha256': model_hash, 'source_hashes': {f'{d}.json': sha(out/f'{d}.json') for d in ['10028','10049','10076']}}, indent=2)+'\n')


if __name__ == '__main__':
    main()
