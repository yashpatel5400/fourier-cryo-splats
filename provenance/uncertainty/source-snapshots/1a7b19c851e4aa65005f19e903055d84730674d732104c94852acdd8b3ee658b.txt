#!/usr/bin/env python3
"""Known-map pose information; not experimental pose-confidence calibration."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from fourier_splats.uq_cell_moments import cell_fourier_moments
from fourier_splats.uq_data import particle_geometry,VoxelReference
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
PROTOCOL=ROOT/'research/uncertainty/ORACLE-POSE-INFORMATION-PROTOCOL.md'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def cell_pose_jacobian(k,q,ctf,coefficients,box,noise,field_A):
    k,q,ctf=map(np.asarray,(k,q,ctf))
    if k.shape[:-1]!=ctf.shape or q.shape!=(*ctf.shape,2) or noise<=0 or field_A<=0:
        raise ValueError('Matching geometry and positive physical scales required')
    moment=cell_fourier_moments(k,coefficients,box,nthreads=1).conj()
    directions=np.stack([np.cross(k,np.eye(3)[j]) for j in range(3)],axis=-2)
    rotation=-2j*np.pi*np.einsum('nqad,nqd->nqa',directions,moment[...,1:4])*np.pi/180.
    translation=-2j*np.pi*q*moment[...,0,None]/field_A
    jacobian=np.concatenate([rotation,translation],axis=-1)*ctf[...,None]/noise
    return np.concatenate([jacobian.real,jacobian.imag],axis=1)


def local_information(jacobian,relative_rank_tolerance=1e-10):
    rows=[]
    for j in np.asarray(jacobian,float):
        information=j.T@j;eigen=np.linalg.eigvalsh(information)
        threshold=relative_rank_tolerance*max(eigen[-1],np.finfo(float).tiny)
        rank=int(np.sum(eigen>threshold))
        record={'information':information.tolist(),'eigenvalues':eigen.tolist(),'rank':rank,
            'relative_rank_tolerance':relative_rank_tolerance}
        if rank==5:
            covariance=np.linalg.inv(information); covariance=.5*(covariance+covariance.T)
            known_shift=np.linalg.inv(information[:3,:3])
            record.update(covariance=covariance.tolist(),coordinate_sd=np.sqrt(np.diag(covariance)).tolist(),
                rotation_rms_degrees=float(np.sqrt(np.trace(covariance[:3,:3]))),
                shift_rms_A=float(np.sqrt(np.trace(covariance[3:,3:]))),
                rotation_max_direction_sd_degrees=float(np.sqrt(np.linalg.eigvalsh(covariance[:3,:3])[-1])),
                shift_max_direction_sd_A=float(np.sqrt(np.linalg.eigvalsh(covariance[3:,3:])[-1])),
                known_shift_rotation_rms_degrees=float(np.sqrt(np.trace(known_shift))))
        else:
            record.update(covariance=None,reason='Rank deficient: no finite five-parameter covariance reported')
        rows.append(record)
    return rows


def main():
    from audit_uq_grid_refinement import model
    files=[Path(__file__),PROTOCOL,ROOT/'tests/test_oracle_pose_information.py',ROOT/'scripts/audit_uq_grid_refinement.py']
    for p in files:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise ValueError('Commit protocol/source before the diagnostic')
    out=BASE/'oracle-pose-information-v1'
    if out.exists():raise RuntimeError('Preserve all earlier outcomes')
    out.mkdir();summary=[]
    for ds,emd in [('10028','2660'),('10049','6487'),('10076','8434')]:
        path=out/f'{ds}.json';begin=time.perf_counter()
        result={'complete':False,'dataset':ds,'scope':'Known-map, supplied-Gaussian-noise local information; not experimental pose calibration or a nonlinear confidence radius.',
            'units':['degree','degree','degree','angstrom','angstrom'],
            'source_snapshot':source_snapshot(ROOT,Path(__file__),[str(p.relative_to(ROOT)) for p in files[1:]]),
            'signals':[]}
        def save():path.write_text(json.dumps(result,indent=2)+'\n')
        save()
        try:
            fitpath=BASE/f'pilot-selected-fixed-v1/{ds}-pilot_region_1.json'
            fit=json.loads(fitpath.read_text());cfg=fit['config'];target=fit['targets'][0]
            wp=fitpath.with_name(f'{ds}-pilot_region_1-{target["width_fraction_field"]}-weights.npz')
            if sha(wp)!=target['source_weights_sha256']:raise ValueError('Archived weights/noise changed')
            saved=np.load(wp);noise=float(saved['noise_std'])
            g=particle_geometry(ROOT,ds,'inference_half0',radius=12,count=128,seed=cfg['seed'])
            np.testing.assert_array_equal(g['indices'],saved['indices'])
            refpath=ROOT/f'data/uncertainty/references/emd_{emd}.map'
            reference=VoxelReference.from_mrc(refpath,box=64).volume.ravel();reference/=np.linalg.norm(reference)
            pilotpath=BASE/f'representation/{ds}/real_particles-spacing-2.0.npz'
            op,pilot,_,check_noise=model(g,np.load(pilotpath),24,noise=noise);pilot=op.expand(pilot)
            np.testing.assert_allclose(np.linalg.norm(pilot),1.,rtol=1e-12);assert noise==check_noise
            result.update(indices=g['indices'].tolist(),field_A=g['field_A'],noise_std_supplied=noise,
                input_hashes={str(p.relative_to(ROOT)):sha(p) for p in [fitpath,wp,refpath,pilotpath,ROOT/f'data/{ds}/metadata.npz',ROOT/f'research/uncertainty/splits/{ds}.csv']})
            for name,signal,box in [('reference_generator',reference,64),('gaussian_pilot_cells',pilot,24)]:
                jacobian=cell_pose_jacobian(g['k'],g['q'],g['ctf'],signal,box,noise,g['field_A'])
                rows=local_information(jacobian)
                full=[r for r in rows if r['rank']==5]
                fields=['rotation_rms_degrees','shift_rms_A','rotation_max_direction_sd_degrees',
                        'shift_max_direction_sd_A','known_shift_rotation_rms_degrees']
                quantiles={field:np.quantile([r[field] for r in full],[0,.1,.5,.9,1.]).tolist() if full else None for field in fields}
                record={'signal':name,'box':box,'unit_continuous_l2_norm':float(np.linalg.norm(signal)),
                    'full_rank_particles':len(full),'particles':rows,'quantile_probabilities':[0,.1,.5,.9,1.],
                    'quantiles':quantiles,'jacobian_sha256':hashlib.sha256(np.asarray(jacobian,dtype='<f8').tobytes()).hexdigest()}
                result['signals'].append(record);save()
                summary.append({'dataset':ds,'signal':name,'full_rank_particles':len(full),'quantiles':quantiles})
                print(ds,name,len(full),quantiles,flush=True)
            result.update(complete=True,seconds=time.perf_counter()-begin);save()
        except Exception as error:
            result.update(complete=True,error=repr(error),seconds=time.perf_counter()-begin);save();raise
    (out/'summary.json').write_text(json.dumps({'complete':True,'records':summary,
        'scope':'Six known-map information diagnostics. These are not experimental pose confidence radii.',
        'source_hashes':{p.name:sha(p) for p in out.glob('*.json')}},indent=2)+'\n')


if __name__=='__main__':main()
