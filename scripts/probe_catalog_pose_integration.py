#!/usr/bin/env python3
"""One fixed-budget catalogue repair of the failed orientation-integral gate."""
import argparse
import json
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.special import softmax
from scipy.spatial.transform import Rotation
from fourier_splats.uq_pose_importance import PoseProposal, importance_summary
from fourier_splats.uq_pose_catalog import CatalogPoseProposal
from probe_adaptive_pose_integration import ROOT, BASE, sha, training_rotations, physical_means
from probe_uq_bispectrum_poses import direct_cells

SOURCES = ['scripts/probe_catalog_pose_integration.py',
           'scripts/probe_adaptive_pose_integration.py',
           'src/fourier_splats/uq_pose_importance.py',
           'src/fourier_splats/uq_pose_catalog.py',
           'research/uncertainty/CATALOG-POSE-INTEGRATION-PROTOCOL.md']


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--dataset', required=True, choices=['10028','10049','10076'])
    ds=parser.parse_args().dataset
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit source and protocol before outcomes')
    out=BASE/'catalog-pose-integration-v1'/ds
    if out.exists():
        raise ValueError('Preserve all numerical gate attempts')
    out.mkdir(parents=True);(out/'images').mkdir();start=time.perf_counter()
    first=BASE/'adaptive-pose-integration-v1'/ds/'summary.json'
    original=json.loads(first.read_text());assert original['complete'] and len(original['cases'])==128
    prior=BASE/'candidate-fisher-score-v1'/ds
    prior_summary=json.loads((prior/'summary.json').read_text())
    assert sha(prior/'arrays.npz')==prior_summary['array']['sha256']
    with np.load(prior/'arrays.npz') as f:
        q=f['q'];transfer=f['transfer'];rho=f['candidate_density'];region=f['region_density']
        training=f['training_candidate_means'][:32768]-.125*f['training_region_means'][:32768]
        first_rotation=f['stage0_first_rotation']
    catalog=training_rotations(ds)
    np.testing.assert_allclose(Rotation.from_quat(catalog[0]).as_matrix(),first_rotation,atol=1e-14)
    plane=np.pad(q,((0,0),(0,1)))
    train_real=np.ascontiguousarray(np.concatenate([training.real,training.imag],axis=1))
    train_norm=np.sum(abs(training)**2,axis=1)
    result=dict(complete=False,dataset=ds,images=128,cases=[],quadrature_levels=[4096,8192],
        source_hashes={s:sha(ROOT/s) for s in SOURCES},git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [first,prior/'summary.json',prior/'arrays.npz']},
        catalogue_temperature=2.,kernel_rotation_std_degrees=8.,component_masses=[.05,.45,.5],
        local_modes_reused_without_refitting=True)
    catalog_path=out/'catalogue-quaternions.npy';np.save(catalog_path,catalog)
    result['catalogue']=dict(path=str(catalog_path.relative_to(ROOT)),sha256=sha(catalog_path),bytes=catalog_path.stat().st_size)
    def save():(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    save()
    for case in original['cases']:
        position=case['position'];state=case['state'];tick=time.perf_counter()
        source=ROOT/case['arrays']['path'];assert sha(source)==case['arrays']['sha256']
        with np.load(source) as f:y=f['observation']
        saved=case['proposal']
        local=PoseProposal(saved['centers'],saved['covariances'],saved['weights'])
        scores=train_real@np.r_[y.real,y.imag]-.5*train_norm
        weights=softmax(scores/2.)
        proposal=CatalogPoseProposal(local,catalog,weights,np.deg2rad(8.))
        row=dict(position=position,index=case['index'],state=state,source_image_sha256=sha(source),
            banks=[],proposal_construction_seconds=time.perf_counter()-tick,
            catalogue_weight_ess=float(1/np.sum(weights**2)))
        allq=[];allp=[];allk=[];allfamily=[];allcomponent=[]
        for bank in [0,1]:
            seed=261020+int(ds)+position*100+state*10+bank
            qq,family,component=proposal.sample(8192,np.random.default_rng(seed))
            tick=time.perf_counter();logp=proposal.log_density(qq);density_seconds=time.perf_counter()-tick
            kernels=np.empty((8192,2));tick=time.perf_counter()
            for begin in range(0,8192,1024):
                end=begin+1024;means=physical_means(qq[begin:end],plane,transfer,rho,region)
                for model in [0,1]:
                    kernels[begin:end,model]=-.5*np.sum(abs(y-(means[0]-.25*model*means[1]))**2,axis=1)
                if position==0 and state==0 and begin==0:
                    k=plane@Rotation.from_quat(qq[0]).as_matrix()
                    direct=np.array([direct_cells(d.ravel(),k)*transfer for d in [rho,region]])
                    error=float(abs(direct-means[:,0]).max())
                    if error>1e-8:raise ArithmeticError('Physical direct-cell check failed')
                    row.setdefault('direct_mean_checks',[]).append(dict(bank=bank,maximum_error=error))
            record=dict(bank=bank,seed=seed,density_seconds=density_seconds,integration_seconds=time.perf_counter()-tick,
                levels=[dict(draws=n,**importance_summary(kernels[:n],logp[:n])) for n in [4096,8192]])
            row['banks'].append(record);allq.append(qq);allp.append(logp);allk.append(kernels)
            allfamily.append(family);allcomponent.append(component)
        path=out/'images'/f'{position:03d}-{state}.npz'
        np.savez_compressed(path,quaternions=np.array(allq),log_proposal=np.array(allp),
            log_kernels=np.array(allk),family=np.array(allfamily),component=np.array(allcomponent),
            catalogue_weights=weights,observation=y)
        row['arrays']=dict(path=str(path.relative_to(ROOT)),sha256=sha(path),bytes=path.stat().st_size)
        row['between_bank_logratio_difference']=abs(row['banks'][0]['levels'][-1]['log_ratio']-row['banks'][1]['levels'][-1]['log_ratio'])
        result['cases'].append(row);save()
        print(ds,position,state,'difference',round(row['between_bank_logratio_difference'],6),'seconds',round(time.perf_counter()-start,2),flush=True)
    differences=np.array([r['between_bank_logratio_difference'] for r in result['cases']])
    medians=np.array([[np.median([r['banks'][b]['levels'][-1]['ess'][m] for r in result['cases']]) for m in [0,1]] for b in [0,1]])
    result.update(complete=True,seconds=time.perf_counter()-start,gate=dict(
        fraction_ratio_difference_at_most_001=float(np.mean(differences<=.01)),median_ess=medians.tolist(),
        passed=bool(np.mean(differences<=.01)>=.9 and np.all(medians>=256))))
    save();print('COMPLETE',ds,result['gate'],flush=True)


if __name__=='__main__':
    main()
