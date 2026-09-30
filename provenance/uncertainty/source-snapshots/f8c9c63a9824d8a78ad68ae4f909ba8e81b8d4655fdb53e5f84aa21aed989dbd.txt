#!/usr/bin/env python3
"""Prospective continuous-orientation computation, not density calibration."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.special import logsumexp
from fourier_splats.uq_continuous_mixture import FourierGaussianOrbit, continuous_mixture_upper
from fourier_splats.uq_physics import evaluate_pairs, gaussian_pair_gram
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
PROTOCOL=ROOT/'research/uncertainty/CONTINUOUS-MIXTURE-PROTOCOL.md'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    files=[Path(__file__),PROTOCOL,
        ROOT/'research/uncertainty/CONTINUOUS-MIXTURE-THEORY.md',
        ROOT/'src/fourier_splats/uq_continuous_mixture.py',
        ROOT/'src/fourier_splats/uq_mixture_validation.py',
        ROOT/'src/fourier_splats/uq_physics.py',ROOT/'src/fourier_splats/uq_data.py',
        ROOT/'tests/test_uq_continuous_mixture.py']
    for p in files:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise ValueError('Commit protocol and implementation before outcomes')
    out=BASE/'continuous-mixture-v1'
    if out.exists():raise RuntimeError('Preserve earlier outcomes')
    out.mkdir();summaries=[];start=time.perf_counter()
    for ds in ['10028','10049','10076']:
        path=out/f'{ds}.json';started=time.perf_counter()
        record=dict(complete=False,dataset=ds,seed=910001+int(ds),
            scope='Full SO(3), Gaussian-pilot generator, known Gaussian noise/CTFs and zero shifts; likelihood computation only, no rejection or experimental coverage claim.',
            source_snapshot=source_snapshot(ROOT,Path(__file__),[str(p.relative_to(ROOT)) for p in files[1:]]),progress=[])
        def save():path.write_text(json.dumps(record,indent=2)+'\n')
        save()
        try:
            info_path=BASE/f'oracle-pose-information-v1/{ds}.json'
            info=json.loads(info_path.read_text())
            fit_path=BASE/f'pilot-selected-fixed-v1/{ds}-pilot_region_1.json'
            cfg=json.loads(fit_path.read_text())['config']
            g=particle_geometry(ROOT,ds,'inference_half0',radius=12,count=128,seed=cfg['seed'])
            np.testing.assert_array_equal(g['indices'],info['indices'])
            ck_path=BASE/f'representation/{ds}/real_particles-spacing-2.0.npz'
            ck=np.load(ck_path);coefficient=ck['coefficients'].copy()
            original_norm=float(np.sqrt(coefficient@gaussian_pair_gram(ck['centers'],ck['sigma'])@coefficient))
            coefficient/=original_norm
            plane=np.pad(g['q'][0],((0,0),(0,1)))
            transfer=g['ctf']/info['noise_std_supplied']
            rng=np.random.default_rng(record['seed'])
            rotations=Rotation.random(128,random_state=rng)
            k=np.einsum('qi,nij->nqj',plane,rotations.as_matrix())
            signal=evaluate_pairs(k,ck['centers'],ck['sigma'],coefficient)*transfer
            signal=np.concatenate([signal.real,signal.imag],axis=1)
            observed=signal+rng.normal(size=signal.shape)
            record.update(particles=128,real_dimension=observed.shape[1],frequency_radius=12,
                field_A=g['field_A'],supplied_noise_std=info['noise_std_supplied'],
                pilot_original_continuous_l2_norm=original_norm,
                mean_squared_whitened_signal=float(np.mean(np.sum(signal**2,axis=1))),
                indices=g['indices'].tolist(),
                input_hashes={str(p.relative_to(ROOT)):sha(p) for p in [info_path,fit_path,ck_path,
                    ROOT/f'data/{ds}/metadata.npz',ROOT/f'research/uncertainty/splits/{ds}.csv']})
            save()
            orbit=FourierGaussianOrbit(plane,ck['centers'],ck['sigma'],coefficient,transfer,observed)
            def progress(row):
                record['progress'].append(row);save();print(ds,row,flush=True)
            result=continuous_mixture_upper(orbit,initial_bins=(4,2,4),max_splits=8192,
                refit_every=256,mixture_iterations=100,tolerance=1.,wall_seconds=900.,callback=progress)
            oracle_angles=rotations.as_euler('ZYZ')
            oracle_kernels=np.array([orbit.log_kernel(a) for a in oracle_angles])
            oracle_lower=float(np.sum(logsumexp(oracle_kernels,axis=0)-np.log(128))+orbit.normalization)
            if oracle_lower>result['log_likelihood_upper']+1e-7:
                raise ArithmeticError('Continuous bound excludes a feasible oracle-support likelihood')
            array_path=out/f'{ds}-arrays.npz'
            arrays={key:value for key,value in result.items() if isinstance(value,np.ndarray)}
            arrays.update(plane=plane,centers=ck['centers'],sigma=ck['sigma'],coefficients=coefficient,
                transfer=transfer,observed=observed,signal=signal,rotations=rotations.as_matrix(),
                oracle_angles=oracle_angles,oracle_log_kernels=oracle_kernels)
            np.savez_compressed(array_path,**arrays)
            scalar={key:value for key,value in result.items() if not isinstance(value,np.ndarray)}
            record.update(complete=True,scientific_run_complete=True,fit=scalar,
                oracle_feasible_log_likelihood=oracle_lower,
                gap_above_best_recorded_feasible=float(result['log_likelihood_upper']-max(oracle_lower,result['log_likelihood_lower'])),
                arrays_file=str(array_path.relative_to(ROOT)),arrays_sha256=sha(array_path),seconds=time.perf_counter()-started)
            save();summaries.append(dict(dataset=ds,fit=scalar,
                oracle_feasible_log_likelihood=oracle_lower,result_sha256=sha(path)))
        except Exception as error:
            record.update(complete=True,scientific_run_complete=False,error=repr(error),seconds=time.perf_counter()-started);save();raise
    (out/'summary.json').write_text(json.dumps(dict(complete=True,records=summaries,
        seconds=time.perf_counter()-start,scope='Continuous denominator computation; no predictive or uncertainty-calibration outcome.'),indent=2)+'\n')


if __name__=='__main__':main()
