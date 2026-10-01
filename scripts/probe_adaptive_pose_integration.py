#!/usr/bin/env python3
"""Frozen classical importance-sampling gate; this does not construct a UQ test."""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from scipy.special import softmax
from fourier_splats.uq_continuous import CellObservationOperator
from fourier_splats.uq_pose_cone_probe import InterpolatedCellFourier
from fourier_splats.uq_pose_importance import PoseProposal,importance_summary
from probe_uq_bispectrum_poses import direct_cells

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_adaptive_pose_integration.py','src/fourier_splats/uq_pose_importance.py',
    'src/fourier_splats/uq_pose_cone_probe.py','research/uncertainty/ADAPTIVE-POSE-INTEGRATION-PROTOCOL.md']


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def training_rotations(ds):
    rng=np.random.default_rng(261017+int(ds));rows=[]
    for _ in range(64):
        rows.append(Rotation.random(512,random_state=rng).as_quat())
        rng.normal(size=(512,220));rng.normal(size=(512,220))
    return np.concatenate(rows)


def start_indices(scores,quaternions):
    order=np.argsort(scores)[::-1];chosen=[];cutoff=np.cos(np.deg2rad(10)/2)
    for index in order:
        if all(abs(quaternions[index]@quaternions[j])<cutoff for j in chosen):chosen.append(int(index))
        if len(chosen)==8:return chosen
    raise ArithmeticError('Fewer than eight separated starts')


def finite_hessian(f,step=.001):
    eye=np.eye(3)*step;zero=np.zeros(3);base=f(zero);h=np.empty((3,3))
    for i in range(3):
        h[i,i]=(f(eye[i])-2*base+f(-eye[i]))/step**2
        for j in range(i):
            h[i,j]=h[j,i]=(f(eye[i]+eye[j])-f(eye[i]-eye[j])-f(-eye[i]+eye[j])+f(-eye[i]-eye[j]))/(4*step**2)
    return h


def physical_means(quaternions,plane,transfer,rho,region):
    rotations=Rotation.from_quat(quaternions).as_matrix()
    k=np.einsum('qi,nij->nqj',plane,rotations);n,nq=k.shape[:2]
    op=CellObservationOperator(k,np.ones((n,nq)),64,1.)
    values=[]
    for density in [rho,region]:
        z=op.forward(density.ravel()).reshape(n,2*nq)
        values.append((z[:,:nq]+1j*z[:,nq:])*transfer)
    return np.array(values)


def fit_proposal(y,guide,plane,transfer,starts,bank_quaternions):
    centers=[];covariances=[];records=[];logmass=[]
    for index in starts:
        base=Rotation.from_quat(bank_quaternions[index])
        def loss(v,base=base):
            k=plane@(base*Rotation.from_rotvec(v)).as_matrix()
            m=guide.values(k)*transfer
            return .5*float(np.sum(abs(y-m)**2))
        first=loss(np.zeros(3))
        fit=minimize(loss,np.zeros(3),method='BFGS',jac='3-point',options=dict(maxiter=120,gtol=1e-5))
        accepted=bool(np.isfinite(fit.fun) and fit.fun<=first)
        center=base*Rotation.from_rotvec(fit.x if accepted else np.zeros(3))
        def local_loss(v):
            return .5*float(np.sum(abs(y-guide.values(plane@(center*Rotation.from_rotvec(v)).as_matrix())*transfer)**2))
        h=finite_hessian(local_loss);e,u=np.linalg.eigh(h)
        precision=np.clip(e,1/np.deg2rad(30)**2,1/np.deg2rad(.25)**2)
        covariance=(u/precision)@u.T
        value=local_loss(np.zeros(3));centers.append(center.as_quat());covariances.append(covariance)
        logmass.append(-value+.5*np.linalg.slogdet(covariance)[1])
        records.append(dict(start_index=index,initial_loss=first,terminal_loss=float(fit.fun),used_loss=value,
            accepted_terminal=accepted,success=bool(fit.success),status=int(fit.status),message=str(fit.message),
            iterations=int(fit.nit),evaluations=int(fit.nfev),hessian_eigenvalues=e.tolist(),
            clipped_rotation_std_degrees=np.rad2deg(1/np.sqrt(precision)).tolist()))
    weights=softmax(logmass)
    # Avoid numerical zero mixture weights while retaining all eight starts.
    weights=np.maximum(weights,np.finfo(float).tiny);weights/=weights.sum()
    proposal=PoseProposal(centers,covariances,weights)
    return proposal,dict(centers=np.array(centers).tolist(),covariances=np.array(covariances).tolist(),weights=weights.tolist(),modes=records)


def main():
    p=argparse.ArgumentParser();p.add_argument('--dataset',required=True,choices=['10028','10049','10076']);ds=p.parse_args().dataset
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():raise ValueError('Commit protocol and source before outcomes')
    out=BASE/'adaptive-pose-integration-v1'/ds
    if out.exists():raise ValueError('Preserve all numerical-gate attempts')
    out.mkdir(parents=True);(out/'images').mkdir();started=time.perf_counter()
    stage=BASE/'matched-information-ledger-v1'/ds;prior=BASE/'candidate-fisher-score-v1'/ds
    a=json.loads((stage/'summary.json').read_text());j=json.loads((prior/'summary.json').read_text())
    assert a['complete'] and j['complete']
    assert sha(stage/'replayed-images.npz')==a['replay']['sha256']
    assert sha(prior/'arrays.npz')==j['array']['sha256']
    result=dict(complete=False,dataset=ds,images=128,quadrature_levels=[4096,8192],cases=[],
        source_hashes={s:sha(ROOT/s) for s in SOURCES},git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [stage/'summary.json',stage/'replayed-images.npz',prior/'summary.json',prior/'arrays.npz']})
    def save():(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    save()
    with np.load(stage/'replayed-images.npz') as f:means=f['means'];noise=f['noise']
    with np.load(prior/'arrays.npz') as f:
        q=f['q'];transfer=f['transfer'];rho=f['candidate_density'];region=f['region_density']
        train=f['training_candidate_means'][:32768]-.125*f['training_region_means'][:32768]
        first_rotation=f['stage0_first_rotation']
    quats=training_rotations(ds)
    np.testing.assert_allclose(Rotation.from_quat(quats[0]).as_matrix(),first_rotation,atol=1e-14)
    plane=np.pad(q,((0,0),(0,1)));guide=InterpolatedCellFourier(rho-.125*region,64,256)
    train_real=np.ascontiguousarray(np.concatenate([train.real,train.imag],axis=1));train_norm=np.sum(abs(train)**2,axis=1)
    indices=np.linspace(0,8191,64,dtype=int);result['source_image_indices']=indices.tolist();save()
    for position,index in enumerate(indices):
        for state in [0,1]:
            y=means[0,index]-.25*state*means[1,index]+noise[index]
            scores=train_real@np.r_[y.real,y.imag]-.5*train_norm
            starts=start_indices(scores,quats);tick=time.perf_counter()
            proposal,fit_record=fit_proposal(y,guide,plane,transfer,starts,quats)
            row=dict(position=position,index=int(index),state=state,proposal=fit_record,fit_seconds=time.perf_counter()-tick,banks=[])
            allq=[];allp=[];allk=[];alllabels=[]
            for bank in [0,1]:
                rng=np.random.default_rng(261020+int(ds)+position*100+state*10+bank)
                qq,labels=proposal.sample(8192,rng);logp=proposal.log_density(qq);kernels=np.empty((8192,2));tick=time.perf_counter()
                for begin in range(0,8192,1024):
                    end=begin+1024;mm=physical_means(qq[begin:end],plane,transfer,rho,region)
                    for model in [0,1]:kernels[begin:end,model]=-.5*np.sum(abs(y-(mm[0]-.25*model*mm[1]))**2,axis=1)
                    if position==0 and state==0 and begin==0:
                        k=plane@Rotation.from_quat(qq[0]).as_matrix()
                        direct=np.array([direct_cells(d.ravel(),k)*transfer for d in [rho,region]])
                        error=float(np.max(abs(direct-mm[:,0])))
                        if error>1e-8:raise ArithmeticError('Independent direct physical means disagree')
                        row.setdefault('direct_mean_checks',[]).append(dict(bank=bank,maximum_error=error))
                record=dict(bank=bank,seed=261020+int(ds)+position*100+state*10+bank,integration_seconds=time.perf_counter()-tick,levels=[])
                for level in [4096,8192]:record['levels'].append(dict(draws=level,**importance_summary(kernels[:level],logp[:level])))
                row['banks'].append(record);allq.append(qq);allp.append(logp);allk.append(kernels);alllabels.append(labels)
            path=out/'images'/f'{position:03d}-{state}.npz'
            np.savez_compressed(path,quaternions=np.array(allq),log_proposal=np.array(allp),log_kernels=np.array(allk),component=np.array(alllabels),observation=y)
            row['arrays']=dict(path=str(path.relative_to(ROOT)),sha256=sha(path),bytes=path.stat().st_size)
            row['between_bank_logratio_difference']=abs(row['banks'][0]['levels'][-1]['log_ratio']-row['banks'][1]['levels'][-1]['log_ratio'])
            result['cases'].append(row);save()
            print(ds,'image',position,'state',state,'difference',round(row['between_bank_logratio_difference'],6),'seconds',round(time.perf_counter()-started,2),flush=True)
    difference=np.array([r['between_bank_logratio_difference'] for r in result['cases']])
    medians=np.array([[np.median([r['banks'][b]['levels'][-1]['ess'][m] for r in result['cases']]) for m in [0,1]] for b in [0,1]])
    result.update(complete=True,seconds=time.perf_counter()-started,gate=dict(fraction_ratio_difference_at_most_001=float(np.mean(difference<=.01)),
        median_ess=medians.tolist(),passed=bool(np.mean(difference<=.01)>=.9 and np.all(medians>=256))),
        optimizer_successes=sum(m['success'] for r in result['cases'] for m in r['proposal']['modes']),optimizer_attempts=1024)
    save();print('COMPLETE',ds,result['seconds'],result['gate'],flush=True)


if __name__=='__main__':main()
