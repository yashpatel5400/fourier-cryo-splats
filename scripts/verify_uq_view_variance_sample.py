#!/usr/bin/env python3
"""Verifier-selected 1% replay using separable exact cell sums and independent moments."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_data import VoxelReference
from probe_uq_bispectrum_poses import direct_cells
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def separable(density,k):
    x=(np.arange(64)+.5)/64-.5
    ex=np.exp(-2j*np.pi*k[:,0,None]*x);ey=np.exp(-2j*np.pi*k[:,1,None]*x);ez=np.exp(-2j*np.pi*k[:,2,None]*x)
    by_x=np.einsum('zyx,qx->qzy',density.reshape(64,64,64),ex,optimize=True)
    by_y=np.einsum('qzy,qy->qz',by_x,ey,optimize=True)
    return np.einsum('qz,qz->q',by_y,ez)*np.prod(np.sinc(k/64),axis=1)/64**1.5


def own_score(y,triads,w):
    n=y.shape[-1];answer=(abs(y)**2/2-1)@w[:n]
    if len(triads):
        a,b,c=triads.T;z=y[:,a]*y[:,b]*y[:,c].conj()/2
        answer+=z.real@w[n:n+len(a)]+z.imag@w[n+len(a):]
    return answer


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'provenance/uncertainty/view-variance-one-percent-verification.json');args=parser.parse_args()
    path=args.output
    if path.exists():raise ValueError('Preserve previous verification')
    start=time.perf_counter();result=dict(complete=False,verifier_seed=202610015,input_hashes={},datasets=[],polynomials=0,group_count_checks=0,maximum_tensor_vs_slow_sum_error=0.)
    h=json.loads((BASE/'bispectrum-hull-v1/summary.json').read_text())
    for ds,emd in [('10028','2660'),('10049','6487')]:
        sp=BASE/f'bispectrum-view-variance-v1/{ds}/summary.json';r=json.loads(sp.read_text());assert r['complete']
        ap=ROOT/r['array']['path'];dp=ROOT/h['arrays'][ds]['path'];mp=ROOT/f'data/uncertainty/references/emd_{emd}.map';gp=BASE/f'mixture-validation-preflight-v2/{ds}.json'
        for p in [sp,ap,dp,mp,gp]:result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        assert sha(ap)==r['array']['sha256'] and sha(dp)==h['arrays'][ds]['sha256']
        with np.load(ap) as f:saved={k:f[k].copy() for k in f.files if k.startswith('cal_')}
        with np.load(dp) as f:d={k:f[k].copy() for k in f.files}
        q=d['q'];plane=np.pad(q,((0,0),(0,1)));transfer=d['transfer'];nq=len(q)
        geom=json.loads(gp.read_text());rho=VoxelReference.from_mrc(mp,box=64).volume.ravel();rho/=np.linalg.norm(rho)
        x=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij');xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
        removed=rho*(1-np.exp(-.5*np.sum((xyz-np.array(geom['target']['center_fraction_field']))**2,axis=1)/(20/geom['field_A'])**2))
        m=r['calibration_views'];indices=np.sort(np.random.default_rng(result['verifier_seed']+int(ds)).choice(m,size=int(np.ceil(.01*m)),replace=False))
        rng=np.random.default_rng(r['stages'][0]['seed']);batch=r['stages'][0]['batch'];selected=set(indices.tolist());checked=0
        edges=np.linspace(.9,1.1,r['amplitude_cells']+1);amps=np.array([-1.,0.,1.,2.]);matrix=np.vander(amps,4,increasing=True)
        for begin in range(0,m,batch):
            n=min(batch,m-begin);rotations=Rotation.random(n,random_state=rng).as_matrix()
            shape=(n,2,r['replicas_per_group'],nq);noise=rng.normal(size=shape)+1j*rng.normal(size=shape)
            for index in sorted(selected.intersection(range(begin,begin+n))):
                local=index-begin;k=plane@rotations[local];e=noise[local].reshape(-1,nq)
                for label,density in [('true',rho),('removed',removed)]:
                    mean=separable(density,k)*transfer
                    if checked==0:
                        slow=direct_cells(density,k[:5])*transfer[:5];err=float(np.max(abs(slow-mean[:5])))
                        if err>1e-8:raise ArithmeticError('Tensor sum mismatch with literal cell sum')
                        result['maximum_tensor_vs_slow_sum_error']=max(result['maximum_tensor_vs_slow_sum_error'],err)
                    for j,c in enumerate(r['cases']):
                        triads=d['triads'] if c['key'].startswith('power_bispectrum') else np.empty((0,3),int);w=d[c['key']+'_direction']
                        values=np.stack([own_score(a*mean+e,triads,w) for a in amps],axis=1)
                        coeff=np.linalg.solve(matrix,values.T).T;events=np.zeros((len(e),len(edges)-1),bool)
                        for ell,coef in enumerate(coeff):
                            p=np.polynomial.Polynomial(coef);roots=p.deriv().roots()
                            for cell,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
                                points=[lo,hi]+[float(z.real) for z in roots if abs(z.imag)<1e-9 and lo<=z.real<=hi]
                                events[ell,cell]=max(p(points))>c['threshold']
                        events=events.reshape(2,r['replicas_per_group'],-1)
                        np.testing.assert_array_equal(events.mean(axis=1).max(axis=1),saved[f'cal_{label}_grouped'][j,index])
                        np.testing.assert_array_equal(events.max(axis=2).mean(axis=1),saved[f'cal_{label}_individual'][j,index])
                        result['polynomials']+=len(e);result['group_count_checks']+=4
                checked+=1
                if checked%64==0:print(ds,'verified views',checked,'seconds',round(time.perf_counter()-start,2),flush=True)
        assert checked==len(indices)
        result['datasets'].append(dict(dataset=ds,views=checked,total_views=m,indices=indices.tolist(),all_selected_group_counts_identical=True))
    result.update(complete=True,seconds=time.perf_counter()-start);path.write_text(json.dumps(result,indent=2)+'\n')
    print('COMPLETE',result['polynomials'],result['group_count_checks'],result['seconds'],flush=True)


if __name__=='__main__':main()
