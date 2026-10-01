#!/usr/bin/env python3
"""Independent first-group polynomial, concentration and critical-value replay."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import binom
from fourier_splats.uq_bispectrum import moment_features

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def radius(x,width,delta):
    average=sum(x)/len(x);squared=sum((z-average)**2 for z in x)/(len(x)-1)
    return np.sqrt(2*squared*np.log(2/delta)/len(x))+7*width*np.log(2/delta)/(3*(len(x)-1))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,
        default=ROOT/'provenance/uncertainty/view-variance-independent-verification.json');args=parser.parse_args()
    if args.output.exists():raise ValueError('Preserve earlier verification')
    result=dict(complete=False,input_hashes={},first_group_checks=[],bound_checks=0,
        event_count_checks=0,critical_value_checks=0,maximum_bound_discrepancy=0.)
    hp=BASE/'bispectrum-hull-v1/summary.json';h=json.loads(hp.read_text())
    oldp=BASE/'bispectrum-monte-carlo-v1/summary.json';old=json.loads(oldp.read_text())
    for p in [hp,oldp]:result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
    for ds in ['10028','10049']:
        sp=BASE/f'bispectrum-view-variance-v1/{ds}/summary.json';r=json.loads(sp.read_text());assert r['complete']
        p=ROOT/r['array']['path'];assert sha(p)==r['array']['sha256']
        ap=ROOT/h['arrays'][ds]['path'];assert sha(ap)==h['arrays'][ds]['sha256']
        for path in [sp,p,ap]:result['input_hashes'][str(path.relative_to(ROOT))]=sha(path)
        with np.load(p) as f:data={k:f[k].copy() for k in f.files}
        with np.load(ap) as f:directions={k:f[k].copy() for k in f.files}
        nq=len(directions['q']);m=r['calibration_views'];replicas=r['replicas_per_group']
        for stage in [0,1]:
            record=r['stages'][stage];rng=np.random.default_rng(record['seed']);batch=record['batch']
            rotations=Rotation.random(batch,random_state=rng).as_matrix()
            np.testing.assert_array_equal(rotations[0],data[f'stage{stage}_example_rotation'])
            shape=(batch,2,replicas,nq) if stage==0 else (batch,nq)
            noise=(rng.normal(size=shape)+1j*rng.normal(size=shape))[0]
            for j,c in enumerate(r['cases']):
                w=directions[c['key']+'_direction']
                triads=directions['triads'] if c['key'].startswith('power_bispectrum') else np.empty((0,3),int)
                for label in ['true','removed']:
                    mean=data[f'stage{stage}_{label}_example_direct'];threshold=c['threshold']
                    if stage==0:
                        e=noise.reshape(-1,nq);amplitudes=np.array([-1.,0.,1.,2.])
                        values=np.stack([moment_features(a*mean+e,triads,1.)@w for a in amplitudes],axis=1)
                        coefficients=np.linalg.solve(np.vander(amplitudes,4,increasing=True),values.T).T
                        events=np.zeros((len(e),r['amplitude_cells']),bool)
                        edges=np.linspace(.9,1.1,r['amplitude_cells']+1)
                        for k,coef in enumerate(coefficients):
                            polynomial=np.polynomial.Polynomial(coef);roots=polynomial.deriv().roots()
                            for cell,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
                                points=[lo,hi]+[float(z.real) for z in roots if abs(z.imag)<1e-9 and lo<=z.real<=hi]
                                events[k,cell]=max(polynomial(points))>threshold
                        events=events.reshape(2,replicas,-1)
                        grouped=events.mean(axis=1).max(axis=1);individual=events.max(axis=2).mean(axis=1)
                        np.testing.assert_array_equal(grouped,data[f'cal_{label}_grouped'][j,0])
                        np.testing.assert_array_equal(individual,data[f'cal_{label}_individual'][j,0])
                        result['first_group_checks'].append(dict(dataset=ds,key=c['key'],map=label,
                            independent_polynomials=len(e),amplitude_cells=r['amplitude_cells'],exact_group_count_match=True))
                    else:
                        score=float(moment_features(mean+noise,triads,1.)@w)
                        np.testing.assert_allclose(score,data[f'heldout_{label}'][j,0],atol=2e-8,rtol=1e-8)
        for j,c in enumerate(r['cases']):
            for label in ['true','removed']:
                x=data[f'cal_{label}_grouped'][j];raw=data[f'cal_{label}_individual'][j]
                assert np.all(x<=raw) and np.all(x>=0) and np.all(raw<=1)
                np.testing.assert_array_equal(x*replicas,np.round(x*replicas))
                np.testing.assert_array_equal(raw*replicas,np.round(raw*replicas))
                oldcase=next(z for z in old['cases'] if z['dataset']==ds and z['key']==c['key'])
                b=oldcase['counts'][f'stage0_{label}_envelope']['frequency']
                saved=c['calibration'][label];assert b==saved['center']
                average=x.mean(axis=1);ordinary=min(1.,average.mean()+radius(average,1.,.001))
                individual=min(1.,raw.mean()+radius(raw.mean(axis=1),1.,.001))
                rad=radius(average,1.,.00025);lo=max(0.,average.mean()-rad);hi=min(1.,average.mean()+rad)
                y=(x[:,0]-b)*(x[:,1]-b);y_width=max(b*b,(1-b)**2)+b*(1-b)
                yu=y.mean()+radius(y,y_width,.0005)
                vu=min(.25,max(0.,yu-max(0.,lo-b,b-hi)**2))
                for row in saved['bounds']:
                    kappa=row['kappa'];expected=[min(1.,kappa*individual),min(1.,kappa*ordinary),
                        min(1.,kappa*hi,hi+np.sqrt((kappa-1)*vu))]
                    actual=[row[n] for n in ['individual_ratio','grouped_ratio','view_variance']]
                    np.testing.assert_allclose(expected,actual,rtol=1e-11,atol=1e-12)
                    result['maximum_bound_discrepancy']=max(result['maximum_bound_discrepancy'],float(max(abs(np.array(actual)-expected))))
                    result['bound_checks']+=3
                scores=data[f'heldout_{label}'][j];count=int(np.sum(scores>c['threshold']))
                assert count==c['heldout'][label]['count'];result['event_count_checks']+=1
            for p in c['projections']:
                n=p['particles'];k=p['reject_at_count'];prob=p['probability_bound']
                assert binom.sf(k-1,n,prob)<=.049+1e-13
                assert k==0 or binom.sf(k-2,n,prob)>.049-1e-13
                result['critical_value_checks']+=1
        print(ds,'independent group polynomials, all bounds/counts/critical values agree',flush=True)
    result['complete']=True;args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('COMPLETE',result['bound_checks'],result['event_count_checks'],result['critical_value_checks'],
        result['maximum_bound_discrepancy'],flush=True)


if __name__=='__main__':main()
