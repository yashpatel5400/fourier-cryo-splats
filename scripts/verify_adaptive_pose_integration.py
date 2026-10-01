#!/usr/bin/env python3
"""Replay saved importance integrals with full 4D covariance and direct cells."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.special import logsumexp
from probe_uq_bispectrum_poses import direct_cells

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def full_covariance_density(q,proposal):
    total=np.full(len(q),.05)
    for center,cov,weight in zip(proposal['centers'],proposal['covariances'],proposal['weights']):
        transform=(Rotation.from_quat(center)*Rotation.from_quat(np.eye(4))).as_quat().T
        sigma=np.eye(4);sigma[:3,:3]=np.asarray(cov)/4;full=transform@sigma@transform.T
        # Solve via its exact orthogonal transform rather than an ill-conditioned global inverse.
        # This representation is independent of the implementation's relative Rotation objects.
        local=q@transform
        quadratic=np.einsum('ni,ij,nj->n',local,np.linalg.inv(sigma),local)
        sign,ld=np.linalg.slogdet(sigma);assert sign==1
        total+=.95*weight*np.exp(-.5*ld)*quadratic**(-2)
    return np.log(total)


def main():
    p=argparse.ArgumentParser();p.add_argument('--dataset',required=True,choices=['10028','10049','10076']);ds=p.parse_args().dataset
    directory=BASE/'adaptive-pose-integration-v1'/ds;sp=directory/'summary.json';j=json.loads(sp.read_text())
    if not j['complete']:raise ValueError('Wait for all declared images, including failures')
    op=directory/'independent-check.json'
    if op.exists():raise ValueError('Preserve every verification attempt')
    result=dict(complete=False,dataset=ds,source_hash=sha(Path(__file__)),summary_hash=sha(sp),images=[],physical_checks=[])
    prior=BASE/'candidate-fisher-score-v1'/ds/'arrays.npz'
    with np.load(prior) as f:q=f['q'];transfer=f['transfer'];rho=f['candidate_density'];region=f['region_density']
    plane=np.pad(q,((0,0),(0,1)));differences=[];esses=[]
    for c in j['cases']:
        path=ROOT/c['arrays']['path'];assert sha(path)==c['arrays']['sha256']
        with np.load(path) as f:quat=f['quaternions'];logp=f['log_proposal'];kernels=f['log_kernels'];y=f['observation']
        assert quat.shape==(2,8192,4) and kernels.shape==(2,8192,2)
        d=full_covariance_density(quat.reshape(-1,4),c['proposal']).reshape(2,8192)
        density_error=float(np.max(abs(d-logp)))
        if density_error>1e-7:raise ArithmeticError('Quaternion density replay failed')
        values=[];ee=[];max_integral_error=0.
        for b in [0,1]:
            for li,n in enumerate([4096,8192]):
                x=kernels[b,:n]-d[b,:n,None];anchor=x.max(0);unnormalized=np.exp(x-anchor)
                integrals=anchor+np.log(unnormalized.mean(0));weights=unnormalized/unnormalized.sum(0)
                ess=1/np.sum(weights*weights,axis=0);lr=float(integrals[1]-integrals[0]);expected=c['banks'][b]['levels'][li]
                err=float(np.max(abs(integrals-expected['log_integrals'])));max_integral_error=max(max_integral_error,err)
                np.testing.assert_allclose(integrals,expected['log_integrals'],atol=1e-8,rtol=0)
                np.testing.assert_allclose(ess,expected['ess'],atol=2e-5,rtol=1e-8)
                np.testing.assert_allclose(weights.max(0),expected['maximum_weight'],atol=1e-8,rtol=0)
                if n==8192:values.append(lr);ee.append(ess)
        differences.append(abs(values[0]-values[1]));esses.append(ee)
        result['images'].append(dict(position=c['position'],state=c['state'],density_error=density_error,log_integral_error=max_integral_error))
        # Physical residual replay at prespecified images, all frequencies, both maps.
        if c['position'] in [0,31,63] and c['state']==0:
            for index in [0,8191]:
                k=plane@Rotation.from_quat(quat[0,index]).as_matrix()
                m=direct_cells(rho.ravel(),k)*transfer;mr=direct_cells(region.ravel(),k)*transfer
                direct=np.array([-.5*np.sum(abs(y-(m-.25*model*mr))**2) for model in [0,1]])
                err=float(np.max(abs(direct-kernels[0,index])))
                if err>1e-7:raise ArithmeticError('Independent direct physical residual failed')
                result['physical_checks'].append(dict(position=c['position'],state=c['state'],bank=0,draw=index,error=err))
        if c['position']%16==0 and c['state']==0:print(ds,'checked',c['position'],flush=True)
    assert len(differences)==128
    medians=np.median(esses,axis=0);fraction=float(np.mean(np.array(differences)<=.01));passed=bool(fraction>=.9 and np.all(medians>=256))
    assert passed==j['gate']['passed'] and fraction==j['gate']['fraction_ratio_difference_at_most_001']
    result.update(complete=True,fraction=fraction,ess_medians=medians.tolist(),passed=passed)
    op.write_text(json.dumps(result,indent=2)+'\n');print('COMPLETE',ds,passed,flush=True)


if __name__=='__main__':main()
