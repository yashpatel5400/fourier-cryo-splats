#!/usr/bin/env python3
"""Declared frame sensitivity of frozen density features, no estimator fitting."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from fourier_splats.uq_data import particle_geometry,particle_observations
from fourier_splats.uq_continuous import cell_forward,cell_target_coefficients
from fourier_splats.uq_continuous_pose import pose_cell_forward
from fourier_splats.uq_group_noise import cell_density_distance
from fourier_splats.uq_intervals import reference_interval_summary
from fourier_splats.uq_pilot_targets import read_target_lock,verify_pilot_scores
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
PROTOCOL=ROOT/'research/uncertainty/REGISTERED-TARGET-SENSITIVITY-PROTOCOL.md'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    for p in [Path(__file__),PROTOCOL]:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise ValueError('Publish protocol and source before outcome calculation')
    out=BASE/'registered-target-sensitivity-v1'
    if out.exists():raise RuntimeError('Retain previous outcomes')
    out.mkdir();summary=[]
    lock_path=ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json'
    lock=read_target_lock(lock_path)
    for locked in lock['datasets']:
        ds=locked['dataset'];start=time.perf_counter();path=out/f'{ds}.json'
        record=dict(complete=False,dataset=ds,conditional=[],experimental=[],cubic=[],
            source_snapshot=source_snapshot(ROOT,Path(__file__),[str(PROTOCOL.relative_to(ROOT)),
                'scripts/audit_uq_grid_refinement.py']),input_hashes={},
            scope='Post-outcome registered-generator sensitivity. Unchanged estimators/intervals; no experimental density calibration.')
        def save():path.write_text(json.dumps(record,indent=2)+'\n')
        def read(p):
            record['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
            r=json.loads(p.read_text())
            if not r.get('complete') or r.get('error'):raise ValueError(f'Unfinished/failed input {p}')
            return r
        def array(p,expected=None):
            digest=sha(p)
            if expected is not None and digest!=expected:raise ValueError(f'Changed array {p}')
            record['input_hashes'][str(p.relative_to(ROOT))]=digest;return np.load(p)
        def budget():
            if time.perf_counter()-start>1800:raise TimeoutError('Declared geometry budget reached')
        save()
        try:
            registration=read(BASE/f'reference-registration-v1/{ds}.json')
            maps=array(ROOT/registration['arrays_file'],registration['arrays_sha256'])
            rho=maps['reference_registered'].ravel().copy();rho/=np.linalg.norm(rho)
            old_rho=maps['reference_original'].ravel().copy();old_rho/=np.linalg.norm(old_rho)
            g=particle_geometry(ROOT,ds,'inference_half0',radius=12,count=128,seed=609315)
            cp=BASE/f'representation/{ds}/real_particles-spacing-2.0.npz'
            ck=array(cp,locked['pilot_checkpoint_sha256'])
            op,pilot,_,_=model(g,ck,24,noise=1.);pilot=op.expand(pilot);verify_pilot_scores(locked,pilot)
            fresh=read(ROOT/f'results/uncertainty/confirmation/noise-calibration-v1/{ds}.json')
            amplitude_record=read(BASE/f'experimental-noise-grouped/{ds}.json')
            gp=particle_geometry(ROOT,ds,'pilot',radius=5,count=256,seed=609681)
            pixels=particle_observations(ROOT,ds,gp);pixels=np.concatenate([pixels.real,pixels.imag],axis=1).ravel()
            pred=cell_forward(gp['k'],gp['ctf'],pilot,24,1.)
            pilot_amplitude=float(pred@pixels/(pred@pred))
            np.testing.assert_allclose(pilot_amplitude,amplitude_record['pilot_amplitude_raw_units'],rtol=1e-10)
            ref_pred=cell_forward(gp['k'],gp['ctf'],rho,64,1.)
            ref_amplitude=float(ref_pred@pixels/(ref_pred@ref_pred))
            if not np.isfinite(ref_amplitude):raise ValueError('Nonfinite reference amplitude')
            scaled_rho=rho*ref_amplitude/pilot_amplitude
            distance,pad=cell_density_distance(pilot,24,scaled_rho,64)
            record['reference_amplitude']=dict(raw=ref_amplitude,pilot_raw=pilot_amplitude,
                ratio=ref_amplitude/pilot_amplitude,positive=ref_amplitude>0,
                reference_pilot_L2_distance=distance,distance_squared_roundoff_pad=pad,
                inside_supplied_radius_two=bool(distance<=2))
            # Noise is the same original supplied scalar for every target.
            first_name=locked['features'][0]['name']
            first=read(BASE/f'pilot-selected-fixed-v1/{ds}-{first_name}.json')
            noise=first['noise_std'];nominal=cell_forward(g['k'],g['ctf'],pilot,24,noise)
            old_signal=cell_forward(g['k'],g['ctf'],old_rho,64,noise)
            signal=cell_forward(g['k'],g['ctf'],rho,64,noise)
            signals={'fixed_pose/nominal':signal}
            rng=np.random.default_rng(610281+int(ds));random=rng.normal(size=(len(g['indices']),5))
            random/=np.linalg.norm(random,axis=1)[:,None]
            for angle in [0,1,2]:
                label='shift_only' if angle==0 else f'pose_{angle}deg'
                for scenario in ['nominal','coherent_x','random_boundary']:
                    budget();poses=np.zeros_like(random)
                    if scenario=='coherent_x':poses[:,0]=1.
                    if scenario=='random_boundary':poses=random.copy()
                    signals[label+'/'+scenario]=pose_cell_forward(g['k'],g['q'],g['ctf'],rho,64,noise,
                        poses,np.deg2rad(angle),.5/g['field_A'])
            replays=[]
            for feature in locked['features']:
                budget();name=feature['name'];centers=[feature['center_fraction_field']];width=locked['width_fraction_field']
                fixed_path=BASE/f'pilot-selected-fixed-v1/{ds}-{name}.json';fixed=read(fixed_path)
                row=fixed['targets'][0];fit=row['fit'];np.testing.assert_allclose(fixed['noise_std'],noise,rtol=1e-14)
                wpath=fixed_path.with_name(f'{ds}-{name}-{width}-weights.npz')
                saved=array(wpath,row['source_weights_sha256']);np.testing.assert_array_equal(saved['indices'],g['indices'])
                w=saved['weights'];pilot_target=feature['pilot_expected_feature'];no_data=2*fit['target_norm']
                old_expected=float(pilot_target+w@(old_signal-nominal))
                np.testing.assert_allclose(old_expected,row['raw_expected_center'],rtol=1e-8,atol=1e-8)
                replays.append(dict(target=name,expected_center_error=old_expected-row['raw_expected_center']))
                target_vector=cell_target_coefficients(64,centers,[1],width);truth=float(target_vector@rho)
                cases=[('fixed_pose',fit['half_width'],['nominal'])]
                for angle in [0,1,2]:
                    e=read(BASE/f'pilot-selected-enclosing-v1/{ds}-{name}-{angle}.json')
                    if e['source_weights_sha256']!=sha(wpath):raise ValueError('Changed enclosing-audit weights')
                    cases.append(('shift_only' if angle==0 else f'pose_{angle}deg',e['half_width'],['nominal','coherent_x','random_boundary']))
                for label,half,scenarios in cases:
                    for scenario in scenarios:
                        expected=float(pilot_target+w@(signals[label+'/'+scenario]-nominal))
                        check=reference_interval_summary(truth,expected,fit['noise_sd'],pilot_target,half,no_data)
                        record['conditional'].append(dict(target=name,pose_class=label,scenario=scenario,
                            raw_expected_center=expected,half_width=half,no_data_half_width=no_data,**check))
                f=next(f for f in fresh['features'] if f['target']==name)
                actual_truth=float(target_vector@scaled_rho)
                for r in f['records']:
                    record['experimental'].append(dict(target=name,pose_class=r['pose_class'],
                        original_reference_target=r['approximate_reference_target'],
                        original_reference_inside=r['approximate_reference_inside_interval'],
                        registered_reference_target=actual_truth,
                        registered_reference_inside=bool(abs(actual_truth-r['interval_center'])<=r['interval_half_width']),
                        interval_center=r['interval_center'],interval_half_width=r['interval_half_width'],
                        excludes_zero=r['excludes_zero'],uses_no_data=r['uses_no_data']))
                if ds=='10049' and name=='pilot_region_1':
                    for family in ['cubic-weight-probe','cubic-coordinate-probe','cubic-subspace-probe']:
                        p=BASE/f'{family}/{ds}-{name}-1.json';c=read(p);a=array(p.with_suffix('.npz'),c['arrays_sha256'])
                        np.testing.assert_array_equal(a['indices'],g['indices']);np.testing.assert_allclose(a['noise_std'],noise)
                        for scenario in ['nominal','coherent_x','random_boundary']:
                            expected=float(c['pilot_target']+a['weights']@(signals['pose_1deg/'+scenario]-nominal))
                            # Saved optimized noise is the affine estimator SD.
                            sd=c['fit']['noise_sd']
                            check=reference_interval_summary(truth,expected,sd,c['pilot_target'],c['refined_half_width'],c['no_data_half_width'])
                            record['cubic'].append(dict(estimator=family,scenario=scenario,raw_expected_center=expected,**check))
                save()
            array_path=out/f'{ds}-arrays.npz'
            np.savez_compressed(array_path,rho_registered_unit=rho,rho_original_unit=old_rho,
                rho_registered_scaled=scaled_rho,indices=g['indices'],k=g['k'],q=g['q'],ctf=g['ctf'],
                nominal=nominal,original_signal=old_signal,random_boundary=random,
                **{key.replace('/','__'):value for key,value in signals.items()})
            record.update(complete=True,scientific_run_complete=True,seconds=time.perf_counter()-start,
                original_nominal_replays=replays,arrays_file=str(array_path.relative_to(ROOT)),arrays_sha256=sha(array_path))
        except Exception as error:
            record.update(complete=True,scientific_run_complete=False,error=repr(error),seconds=time.perf_counter()-start)
        save();summary.append(dict(dataset=ds,scientific_run_complete=record['scientific_run_complete'],
            result_sha256=sha(path),error=record.get('error'),seconds=record['seconds']))
        print(ds,summary[-1],flush=True)
    (out/'summary.json').write_text(json.dumps(dict(complete=True,records=summary),indent=2)+'\n')


if __name__=='__main__':main()
