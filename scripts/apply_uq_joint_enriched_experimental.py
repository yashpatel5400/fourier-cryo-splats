#!/usr/bin/env python3
"""Declared raw/centered application of the enriched and joint estimators."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from fourier_splats.physics import fft_center
from fourier_splats.uq_data import particle_geometry, particle_observations
from fourier_splats.uq_continuous import cell_forward
from fourier_splats.uq_group_noise import grouped_estimator_variance_upper
from fourier_splats.uq_centered_noise import centered_grouped_variance_upper
from fourier_splats.uq_noise_calibration import uncenter_fourier_weights
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_pilot_targets import read_target_lock, verify_pilot_scores
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
PROTOCOL = ROOT/'research/uncertainty/JOINT-ENRICHED-EXPERIMENTAL-PROTOCOL.md'
FAMILIES = ['cubic-enrichment-probe', 'joint-cubic-design-probe']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    for p in [Path(__file__), PROTOCOL]:
        if subprocess.check_output(['git', 'show', f'HEAD:{p.relative_to(ROOT)}'], cwd=ROOT) != p.read_bytes():
            raise ValueError('Commit protocol and source before application')
    out = BASE/'joint-enriched-experimental-v1'
    if out.exists():
        raise RuntimeError('Preserve previous outcomes')
    out.mkdir(); start = time.perf_counter(); path = out/'10049.json'
    record = dict(complete=False, scientific_run_complete=False, dataset='10049', target='pilot_region_1',
        scope='Exploratory reuse of experimental pixels; no experimental coverage claim.',
        source_snapshot=source_snapshot(ROOT, Path(__file__), [str(PROTOCOL.relative_to(ROOT)),
            'scripts/audit_uq_grid_refinement.py']), input_hashes={}, records=[],
        alpha_noise=(.045-1e-6)/12, beta_calibration=.005/12, delta_spectral=1e-6/12,
        calibration_alternatives=['uncentered','centered'],
        multiplicity_scope='Separate post-development procedures; no unadjusted minimum or experimental simultaneous-coverage claim.',
        rotation_radius_degrees=1., translation_radius_A=.5,
        density_radius=2., pilot_norm_bound=1., width_A=20.)

    def save():
        path.write_text(json.dumps(record, indent=2)+'\n')

    def file(p, expected=None):
        digest = sha(p)
        if expected is not None and digest != expected:
            raise ValueError(f'Changed file: {p}')
        record['input_hashes'][str(p.relative_to(ROOT))] = digest
        return p

    def read(p):
        r = json.loads(file(p).read_text())
        if not r.get('complete') or r.get('error'):
            raise ValueError(f'Incomplete input: {p}')
        return r

    save()
    try:
        dp = ROOT/'provenance/uncertainty/noise-calibration-v1/10049-download.json'
        download = read(dp)
        for name, digest in download['cached_files'].items():
            file(ROOT/name, digest)
        selection_path = ROOT/'research/uncertainty/confirmation/noise-calibration-v1/10049-selection.csv'
        selection = list(csv.DictReader(file(selection_path).open()))
        labels = [r['source_group'] for r in selection]
        excluded = set()
        for name in ['research/uncertainty/splits/10049.csv',
                     'research/uncertainty/confirmation/prediction-v1/10049-selection.csv']:
            excluded.update(r['source_group'] for r in csv.DictReader(file(ROOT/name).open()))
        if len(labels) != 128 or len(set(labels)) != 128 or set(labels) & excluded:
            raise ValueError('Exposure disjointness failed')
        g = particle_geometry(ROOT, '10049', 'inference_half0', radius=12, count=128, seed=609315)
        lock_path = file(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
        locked = next(r for r in read_target_lock(lock_path)['datasets'] if r['dataset'] == '10049')
        cp = file(BASE/'representation/10049/real_particles-spacing-2.0.npz', locked['pilot_checkpoint_sha256'])
        op, pilot, _, _ = model(g, np.load(cp), 24, noise=1.)
        pilot = op.expand(pilot); verify_pilot_scores(locked, pilot)
        original = read(BASE/'experimental-noise-grouped/10049.json')
        amplitude = original['pilot_amplitude_raw_units']
        images = np.load(file(ROOT/'data/uncertainty/confirmation/noise-calibration-v1/10049/images.npy'))
        q = g['q'][0].astype(int); c = images.shape[-1]//2
        ft = fft_center(images.astype(float))[:, q[:, 1]+c, q[:, 0]+c]
        calibration = np.concatenate([ft.real, ft.imag], axis=1)/amplitude
        observation = particle_observations(ROOT, '10049', g)
        y = np.concatenate([observation.real, observation.imag], axis=1).ravel()/amplitude
        phases = g['manifest']['data_sign']*np.exp(-2j*np.pi*np.einsum('nqi,ni->nq', g['q'], g['translations']))
        fixed = read(BASE/'pilot-selected-fixed-v1/10049-pilot_region_1.json')
        f = fixed['targets'][0]; width = locked['width_fraction_field']
        wp = file(BASE/f'pilot-selected-fixed-v1/10049-pilot_region_1-{width}-weights.npz', f['source_weights_sha256'])
        saved = np.load(wp); np.testing.assert_array_equal(saved['indices'], g['indices'])
        noise = float(saved['noise_std']); nominal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
        applied = read(BASE/'pilot-selected-experimental-v1/10049-pilot_region_1.json')
        fresh = read(ROOT/'results/uncertainty/confirmation/noise-calibration-v1/10049.json')
        fr = next(r for r in fresh['features'] if r['target'] == 'pilot_region_1')
        fresh_fixed = next(r for r in fr['records'] if r['pose_class'] == 'fixed_pose')
        raw = uncenter_fourier_weights(saved['weights'].reshape(128, -1), phases, noise)
        calibrated = grouped_estimator_variance_upper(calibration, raw, g['groups'], record['beta_calibration'])
        center = float(f['pilot_expected_feature']+saved['weights']@(y/noise-nominal))
        half = bias_aware_half_width_stable(calibrated['noise_sd_upper'], fresh_fixed['bias_upper'], record['alpha_noise'])
        np.testing.assert_allclose([center, calibrated['noise_sd_upper'], half],
            [fresh_fixed['raw_observed_center'], fresh_fixed['noise_sd_upper'], fresh_fixed['raw_half_width']],
            rtol=1e-9, atol=1e-9)
        record['fixed_estimator_replay'] = dict(center=center, noise_sd_upper=calibrated['noise_sd_upper'], half_width=half)
        registered = read(BASE/'registered-target-sensitivity-v1/10049.json')
        rr = next(r for r in registered['experimental'] if r['target']=='pilot_region_1' and r['pose_class']=='pose_1deg')
        record.update(calibration_groups=labels, inference_groups=g['groups'].tolist(),
            calibration_group_count=128, excluded_group_overlap=0, pilot_amplitude_raw_units=amplitude,
            reference_amplitude=registered['reference_amplitude'],
            unverified_conditions=original['unverified_conditions']+[
                calibrated['assumptions'], 'Both inference and calibration data have previously been analyzed.',
                'Pose radii, noise homogeneity, supplied-pose independence, and density radius are unverified.',
                'Reference inclusion measures agreement; the deposited reference is not density truth.'])
        arrays = dict(calibration=calibration, observed=y, nominal=nominal, indices=g['indices'],
            phases=phases, k=g['k'], q=g['q'], ctf=g['ctf'], noise_std=noise)
        # Replay both existing calibration procedures before new outcomes.
        centered_record = read(BASE/'centered-noise-calibration-v1/10049.json')
        old_application = read(BASE/'cubic-experimental-application-v1/10049.json')
        old_arrays = np.load(file(ROOT/old_application['arrays_file'], old_application['arrays_sha256']))
        np.testing.assert_allclose(old_arrays['calibration'],calibration,rtol=1e-12)
        record['old_cubic_calibration_replays'] = []
        for previous in centered_record['cubic']:
            old_raw=old_arrays[previous['estimator']+'_raw_weights']
            raw_cal=grouped_estimator_variance_upper(calibration,old_raw,g['groups'],record['beta_calibration'])
            cen_cal,_,_=centered_grouped_variance_upper(calibration,old_raw,g['groups'],record['beta_calibration'])
            np.testing.assert_allclose([raw_cal['noise_sd_upper'],cen_cal['noise_sd_upper']],
                [previous['calibration']['uncentered']['noise_sd_upper'],previous['calibration']['centered']['noise_sd_upper']],
                rtol=1e-9,atol=1e-9)
            record['old_cubic_calibration_replays'].append(dict(estimator=previous['estimator'],
                uncentered_sd=raw_cal['noise_sd_upper'],centered_sd=cen_cal['noise_sd_upper']))
        save()
        for family in FAMILIES:
            if time.perf_counter()-start > 900:
                raise TimeoutError('Declared 15-minute budget reached')
            source = BASE/f'{family}/10049-pilot_region_1-1.json'; old = read(source)
            a = np.load(file(source.with_suffix('.npz'), old['arrays_sha256']))
            np.testing.assert_array_equal(a['indices'], g['indices'])
            np.testing.assert_allclose(a['noise_std'], noise, rtol=1e-14)
            np.testing.assert_allclose(old['pilot_target'], f['pilot_expected_feature'], rtol=1e-12)
            np.testing.assert_allclose(old['width_fraction_field']*g['field_A'], 20., rtol=1e-12)
            w = a['weights']; raw = uncenter_fourier_weights(w.reshape(128, -1), phases, noise)
            center = float(old['pilot_target']+w@(y/noise-nominal))
            bias = old['bound']['bias_upper'] if family=='cubic-enrichment-probe' else old['fit']['bias_upper']
            np.testing.assert_allclose(bias_aware_half_width_stable(old['fit']['noise_sd'],bias,
                old['fit']['alpha_noise']),old['refined_half_width'],rtol=1e-12)
            raw_cal=grouped_estimator_variance_upper(calibration,raw,g['groups'],record['beta_calibration'])
            centered_cal,contrasts,centered_data=centered_grouped_variance_upper(calibration,raw,g['groups'],record['beta_calibration'])
            arrays[family+'_raw_weights']=raw
            arrays['helmert_contrasts']=contrasts;arrays['centered_calibration']=centered_data
            for method,cal in [('uncentered',raw_cal),('centered',centered_cal)]:
                raw_half=bias_aware_half_width_stable(cal['noise_sd_upper'],bias,record['alpha_noise'])
                no_data=old['no_data_half_width'];fallback=bool(raw_half>=no_data)
                interval_center=old['pilot_target'] if fallback else center;half=min(raw_half,no_data)
                row=dict(estimator=family,calibration_method=method,bias_upper=bias,raw_observed_center=center,
                    variance_calibration=cal,supplied_simulation_noise_sd=old['fit']['noise_sd'],
                    calibrated_to_simulation_sd=cal['noise_sd_upper']/old['fit']['noise_sd'],
                    raw_half_width=raw_half,no_data_half_width=no_data,interval_center=interval_center,
                    interval_half_width=half,relative_half_width=half/no_data,uses_no_data=fallback,
                    excludes_zero=bool(abs(interval_center)>half),
                    original_reference_target=rr['original_reference_target'],registered_reference_target=rr['registered_reference_target'],
                    original_reference_inside=bool(abs(rr['original_reference_target']-interval_center)<=half),
                    registered_reference_inside=bool(abs(rr['registered_reference_target']-interval_center)<=half))
                record['records'].append(row);save()
                print(family,method,'center',center,'half',half,'SD',cal['noise_sd_upper'],flush=True)
        ap = out/'10049-arrays.npz'; np.savez_compressed(ap, **arrays)
        record.update(complete=True, scientific_run_complete=True, seconds=time.perf_counter()-start,
            arrays_file=str(ap.relative_to(ROOT)), arrays_sha256=sha(ap))
    except Exception as error:
        record.update(complete=True, scientific_run_complete=False, error=repr(error), seconds=time.perf_counter()-start)
        save(); raise
    save()


if __name__ == '__main__':
    main()
