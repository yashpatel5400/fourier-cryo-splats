#!/usr/bin/env python3
"""Declared centered-noise sensitivity with all prior estimators fixed."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from fourier_splats.physics import fft_center
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_noise_calibration import uncenter_fourier_weights
from fourier_splats.uq_group_noise import grouped_estimator_variance_upper
from fourier_splats.uq_centered_noise import centered_grouped_variance_upper
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_pilot_targets import read_target_lock
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'
PROTOCOL = ROOT/'research/uncertainty/CENTERED-NOISE-CALIBRATION-PROTOCOL.md'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    for p in [Path(__file__), PROTOCOL, ROOT/'src/fourier_splats/uq_centered_noise.py']:
        if subprocess.check_output(['git', 'show', f'HEAD:{p.relative_to(ROOT)}'], cwd=ROOT) != p.read_bytes():
            raise ValueError('Commit procedure before outcomes')
    out = BASE/'centered-noise-calibration-v1'
    if out.exists(): raise RuntimeError('Preserve previous outcomes')
    out.mkdir(); start = time.perf_counter(); summary = []
    lock = read_target_lock(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
    for locked in lock['datasets']:
        ds = locked['dataset']; path = out/f'{ds}.json'
        r = dict(complete=False, scientific_run_complete=False, dataset=ds, records=[], cubic=[], input_hashes={},
            source_snapshot=source_snapshot(ROOT, Path(__file__), [str(PROTOCOL.relative_to(ROOT))]),
            scope='Reused-data centered-noise sensitivity; not experimental coverage validation.',
            alpha_noise=(.045-1e-6)/12, beta_calibration=.005/12, delta_spectral=1e-6/12)
        def save(): path.write_text(json.dumps(r, indent=2)+'\n')
        def file(p, expected=None):
            digest=sha(p)
            if expected is not None and expected!=digest: raise ValueError(f'Changed input: {p}')
            r['input_hashes'][str(p.relative_to(ROOT))]=digest
            return p
        def read(p):
            data=json.loads(file(p).read_text())
            if not data.get('complete') or data.get('error'): raise ValueError(f'Incomplete source: {p}')
            return data
        def interval(source, sd, pilot_target, no_data):
            half=bias_aware_half_width_stable(sd, source['bias_upper'], r['alpha_noise'])
            fallback=bool(half>=no_data); raw_half=half; half=min(half, no_data)
            center=pilot_target if fallback else source['raw_observed_center']
            return dict(raw_half_width=raw_half, interval_center=center, interval_half_width=half,
                relative_half_width=half/no_data, uses_no_data=fallback, excludes_zero=bool(abs(center)>half),
                original_reference_inside=bool(abs(source['original_reference_target']-center)<=half),
                registered_reference_inside=bool(abs(source['registered_reference_target']-center)<=half))
        def calibration(raw, expected):
            uncentered=grouped_estimator_variance_upper(y, raw, g['groups'], r['beta_calibration'])
            np.testing.assert_allclose(uncentered['noise_sd_upper'], expected, rtol=1e-9, atol=1e-9)
            centered,q,contrasts=centered_grouped_variance_upper(y, raw, g['groups'], r['beta_calibration'])
            energy=uncentered['calibration']['observed_total_energy']; ce=centered['calibration']['observed_total_energy']
            ratio=centered['noise_sd_upper']/expected if expected else None
            report=dict(uncentered=uncentered, centered=centered, centered_to_uncentered_sd=ratio,
                centered_energy_fraction=ce/energy if energy else None)
            arrays['contrast_matrix']=q; arrays['centered_calibration']=contrasts
            return report
        save()
        try:
            file(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
            dp=read(ROOT/f'provenance/uncertainty/noise-calibration-v1/{ds}-download.json')
            for name,digest in dp['cached_files'].items(): file(ROOT/name,digest)
            old=read(BASE/f'experimental-noise-grouped/{ds}.json')
            fresh=read(ROOT/f'results/uncertainty/confirmation/noise-calibration-v1/{ds}.json')
            registered=read(BASE/f'registered-target-sensitivity-v1/{ds}.json')
            g=particle_geometry(ROOT,ds,'inference_half0',radius=12,count=128,seed=609315)
            images=np.load(file(ROOT/f'data/uncertainty/confirmation/noise-calibration-v1/{ds}/images.npy'))
            q=g['q'][0].astype(int); c=images.shape[-1]//2
            ft=fft_center(images.astype(float))[:,q[:,1]+c,q[:,0]+c]
            y=np.concatenate([ft.real,ft.imag],axis=1)/old['pilot_amplitude_raw_units']
            phases=g['manifest']['data_sign']*np.exp(-2j*np.pi*np.einsum('nqi,ni->nq',g['q'],g['translations']))
            arrays=dict(calibration=y,indices=g['indices'],phases=phases)
            r['unverified_conditions']=fresh['unverified_conditions']+[
                'Centered procedure developed after prior calibration outcomes; no fresh confirmation.',
                'Reference agreement is not truth coverage; no minimum over procedures is claimed.']
            for feature in locked['features']:
                if time.perf_counter()-start>900: raise TimeoutError('Declared full-study budget reached')
                name=feature['name']; fp=BASE/f'pilot-selected-fixed-v1/{ds}-{name}.json'; fixed=read(fp)
                f=fixed['targets'][0]; no_data=2*f['fit']['target_norm']
                wp=file(fp.with_name(f'{ds}-{name}-{locked["width_fraction_field"]}-weights.npz'),f['source_weights_sha256'])
                a=np.load(wp); np.testing.assert_array_equal(a['indices'],g['indices'])
                raw=uncenter_fourier_weights(a['weights'].reshape(128,-1),phases,float(a['noise_std']))
                fr=next(row for row in fresh['features'] if row['target']==name)
                cal=calibration(raw,fr['variance_calibration']['noise_sd_upper']); sd=cal['centered']['noise_sd_upper']
                rows=[]
                for previous in fr['records']:
                    file(ROOT/previous['bias_source'],previous['bias_source_sha256'])
                    ref=next(row for row in registered['experimental'] if row['target']==name and row['pose_class']==previous['pose_class'])
                    source=dict(previous,original_reference_target=ref['original_reference_target'],registered_reference_target=ref['registered_reference_target'])
                    rows.append(dict(pose_class=previous['pose_class'],bias_upper=previous['bias_upper'],
                        raw_observed_center=previous['raw_observed_center'],old_interval_half_width=previous['interval_half_width'],
                        old_excludes_zero=previous['excludes_zero'],original_reference_target=ref['original_reference_target'],
                        registered_reference_target=ref['registered_reference_target'],**interval(source,sd,feature['pilot_expected_feature'],no_data)))
                r['records'].append(dict(target=name,calibration=cal,intervals=rows)); arrays[name+'_raw_weights']=raw; save()
            if ds=='10049':
                original=read(BASE/'cubic-experimental-application-v1/10049.json')
                source_arrays=np.load(file(ROOT/original['arrays_file'],original['arrays_sha256']))
                np.testing.assert_allclose(source_arrays['calibration'],y,rtol=1e-12)
                for source in original['records']:
                    family=source['estimator']; raw=source_arrays[family+'_raw_weights']
                    cal=calibration(raw,source['variance_calibration']['noise_sd_upper'])
                    case=read(BASE/f'{family}/10049-pilot_region_1-1.json')
                    r['cubic'].append(dict(estimator=family,calibration=cal,bias_upper=source['bias_upper'],
                        raw_observed_center=source['raw_observed_center'],old_interval_half_width=source['interval_half_width'],
                        original_reference_target=source['original_reference_target'],registered_reference_target=source['registered_reference_target'],
                        **interval(source,cal['centered']['noise_sd_upper'],case['pilot_target'],source['no_data_half_width'])))
                    arrays[family+'_raw_weights']=raw; save()
            ap=out/f'{ds}-arrays.npz'; np.savez_compressed(ap,**arrays)
            r.update(complete=True,scientific_run_complete=True,arrays_file=str(ap.relative_to(ROOT)),arrays_sha256=sha(ap),
                seconds_since_study_start=time.perf_counter()-start)
        except Exception as error:
            r.update(complete=True,scientific_run_complete=False,error=repr(error),seconds_since_study_start=time.perf_counter()-start)
            save(); raise
        save(); summary.append(dict(dataset=ds,sha256=sha(path),scientific_run_complete=True))
        print('DONE',ds,[(row['target'],row['calibration']['centered_to_uncentered_sd']) for row in r['records']],flush=True)
    (out/'summary.json').write_text(json.dumps(dict(complete=True,records=summary),indent=2)+'\n')


if __name__=='__main__': main()
