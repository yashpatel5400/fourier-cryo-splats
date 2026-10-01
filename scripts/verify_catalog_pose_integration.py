#!/usr/bin/env python3
"""All saved integrals; full-4D density checks on 32 fixed draws per bank/image."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from verify_adaptive_pose_integration import ROOT,BASE,sha
from probe_uq_bispectrum_poses import direct_cells


def density_full_covariance(points, local, inv_catalog, determinant_catalog, weights):
    total=np.full(len(points),.05)
    for center,cov,weight in zip(local['centers'],local['covariances'],local['weights']):
        transform=(Rotation.from_quat(center)*Rotation.from_quat(np.eye(4))).as_quat().T
        sigma=np.eye(4);sigma[:3,:3]=np.asarray(cov)/4
        q=points@transform
        z=np.einsum('ni,ij,nj->n',q,np.linalg.inv(sigma),q)
        total+=.45*weight*np.linalg.det(sigma)**(-.5)/z**2
    # Explicit inverse 4x4 matrices; no dot-product ACG formula from the runner.
    z=np.einsum('qi,cij,qj->cq',points,inv_catalog,points,optimize=True)
    values=determinant_catalog[:,None]**(-.5)/z**2
    total+=.5*(weights@values)
    return np.log(total)


def main():
    p=argparse.ArgumentParser();p.add_argument('--dataset',required=True,choices=['10028','10049','10076']);ds=p.parse_args().dataset
    directory=BASE/'catalog-pose-integration-v1'/ds
    j=json.loads((directory/'summary.json').read_text());assert j['complete']
    op=directory/'independent-check.json'
    if op.exists():raise ValueError('Preserve verification attempts')
    previous=json.loads((BASE/'adaptive-pose-integration-v1'/ds/'summary.json').read_text())
    catalog_path=ROOT/j['catalogue']['path'];assert sha(catalog_path)==j['catalogue']['sha256']
    catalog=np.load(catalog_path);s=np.deg2rad(8)/2
    sigma=s*s*np.eye(4)[None]+(1-s*s)*np.einsum('ci,cj->cij',catalog,catalog)
    inverse=np.linalg.inv(sigma);determinant=np.linalg.det(sigma)
    result=dict(complete=False,dataset=ds,source_hash=sha(Path(__file__)),summary_hash=sha(directory/'summary.json'),
        proposal_check_scope='32 integer-linspace draws in each of two banks for every image; not every sampled density',
        checked_draw_indices=np.linspace(0,8191,32,dtype=int).tolist(),images=[],physical_checks=[])
    with np.load(BASE/'candidate-fisher-score-v1'/ds/'arrays.npz') as f:
        q=f['q'];transfer=f['transfer'];rho=f['candidate_density'];region=f['region_density']
    plane=np.pad(q,((0,0),(0,1)));differences=[];esses=[]
    for c,old in zip(j['cases'],previous['cases']):
        assert (c['position'],c['state'],c['index'])==(old['position'],old['state'],old['index'])
        path=ROOT/c['arrays']['path'];assert sha(path)==c['arrays']['sha256']
        with np.load(path) as f:
            quat=f['quaternions'];logp=f['log_proposal'];kernels=f['log_kernels'];y=f['observation'];weights=f['catalogue_weights']
        with np.load(ROOT/old['arrays']['path']) as f:np.testing.assert_array_equal(y,f['observation'])
        np.testing.assert_allclose(np.linalg.norm(quat,axis=-1),1,atol=2e-14,rtol=0)
        ix=result['checked_draw_indices'];points=quat[:,ix].reshape(-1,4)
        independent=density_full_covariance(points,old['proposal'],inverse,determinant,weights).reshape(2,-1)
        density_error=float(abs(independent-logp[:,ix]).max())
        if density_error>1e-7:raise ArithmeticError('Full covariance density check failed')
        ratios=[];ee=[];integral_error=0.
        for bank in [0,1]:
            for li,n in enumerate([4096,8192]):
                x=kernels[bank,:n]-logp[bank,:n,None];anchor=x.max(0)
                u=np.exp(x-anchor);integrals=anchor+np.log(u.mean(0));w=u/u.sum(0)
                ess=1/np.sum(w*w,axis=0);expected=c['banks'][bank]['levels'][li]
                integral_error=max(integral_error,float(abs(integrals-expected['log_integrals']).max()))
                np.testing.assert_allclose(integrals,expected['log_integrals'],atol=1e-10,rtol=0)
                np.testing.assert_allclose(ess,expected['ess'],atol=1e-7,rtol=1e-10)
                np.testing.assert_allclose(w.max(0),expected['maximum_weight'],atol=1e-12,rtol=1e-10)
                ratio=float(integrals[1]-integrals[0]);np.testing.assert_allclose(ratio,expected['log_ratio'],atol=1e-10,rtol=0)
                if n==8192:ratios.append(ratio);ee.append(ess)
        differences.append(abs(ratios[0]-ratios[1]));esses.append(ee)
        np.testing.assert_allclose(differences[-1],c['between_bank_logratio_difference'],atol=1e-10,rtol=0)
        result['images'].append(dict(position=c['position'],state=c['state'],checked_densities=len(points),
            maximum_log_density_error=density_error,maximum_log_integral_error=integral_error))
        if c['position'] in [0,31,63] and c['state']==0:
            for index in [0,8191]:
                k=plane@Rotation.from_quat(quat[0,index]).as_matrix()
                m=direct_cells(rho.ravel(),k)*transfer;mr=direct_cells(region.ravel(),k)*transfer
                direct=np.array([-.5*np.sum(abs(y-(m-.25*model*mr))**2) for model in [0,1]])
                error=float(abs(direct-kernels[0,index]).max())
                if error>1e-7:raise ArithmeticError('Physical cell-sum check failed')
                result['physical_checks'].append(dict(position=c['position'],bank=0,draw=index,error=error))
        if c['position']%16==0 and c['state']==0:print(ds,'checked',c['position'],flush=True)
    assert len(differences)==128
    median=np.median(esses,axis=0);fraction=float(np.mean(np.array(differences)<=.01));passed=bool(fraction>=.9 and np.all(median>=256))
    assert passed==j['gate']['passed'] and fraction==j['gate']['fraction_ratio_difference_at_most_001']
    np.testing.assert_allclose(median,j['gate']['median_ess'],atol=1e-7,rtol=1e-10)
    result.update(complete=True,checked_densities=sum(r['checked_densities'] for r in result['images']),
                  fraction=fraction,ess_medians=median.tolist(),passed=passed)
    op.write_text(json.dumps(result,indent=2)+'\n');print('COMPLETE',ds,passed,flush=True)


if __name__=='__main__':main()
