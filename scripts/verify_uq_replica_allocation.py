#!/usr/bin/env python3
"""Independent inherited-score, new allocation calibration and all-outcome replay."""
import argparse,hashlib,json
from pathlib import Path
import mrcfile
import numpy as np
from scipy.signal import fftconvolve
from scipy.linalg import null_space
from scipy.spatial.transform import Rotation
from scipy.stats import binom,beta
from verify_uq_view_variance_sample import separable,own_score
from verify_uq_view_variance import radius
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def features(m,triads):
    power=(m.real*m.real+m.imag*m.imag)/2
    if len(triads)==0:return power
    products=m[:,triads[:,0]]*m[:,triads[:,1]]*np.conjugate(m[:,triads[:,2]])
    return np.column_stack([power,products.real/2,products.imag/2])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'provenance/uncertainty/replica-allocation-independent-verification.json');args=parser.parse_args()
    if args.output.exists():raise ValueError('Preserve prior verification')
    result=dict(complete=False,input_hashes={},candidate_mask_checks=0,calibration_first_polynomials=0,
        direct_score_checks=0,event_count_checks=0,critical_value_checks=0,group_vector_checks=0,projection_checks=0,
        inherited_score_checks=0,bound_checks=0,maximum_bound_difference=0.,maximum_physical_coordinate_difference=0.,maximum_direct_score_difference=0.)
    for ds in ['10028','10049','10076']:
        sp=BASE/f'candidate-replica-allocation-v1/{ds}/summary.json';r=json.loads(sp.read_text());assert r['complete']
        ap=ROOT/r['array']['path'];assert sha(ap)==r['array']['sha256'];mp=ROOT/f'results/final/{ds}/gaussian-mean.mrc'
        for p in [sp,ap,mp]:result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        with np.load(ap) as f:d={k:f[k].copy() for k in f.files}
        with mrcfile.open(mp) as f:rho=np.array(f.data,float);field=float(f.voxel_size.x)*64
        rho/=np.linalg.norm(rho);np.testing.assert_array_equal(rho,d['candidate_density']);assert field==r['field_A']
        # Independent 3-D FFT convolution against the separable Gaussian filter.
        sigma=20/field*64;half=int(4*sigma+.5);x=np.arange(-half,half+1);kernel=np.exp(-x*x/(2*sigma*sigma));kernel/=kernel.sum()
        k3=kernel[:,None,None]*kernel[None,:,None]*kernel[None,None,:];smooth=fftconvolve(rho,k3,mode='same')
        grid=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(grid,grid,grid,indexing='ij');inside=(abs(x)<=.25)&(abs(y)<=.25)&(abs(z)<=.25)
        peak=np.unravel_index(np.argmax(np.where(inside,smooth,-np.inf)),rho.shape);assert list(peak)==r['region']['center_index_zyx']
        center=np.array([grid[peak[2]],grid[peak[1]],grid[peak[0]]]);mask=np.exp(-.5*((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)/(20/field)**2)
        np.testing.assert_array_equal(mask,d['region_mask']);np.testing.assert_array_equal(rho*mask,d['region_density']);result['candidate_mask_checks']+=1
        q=d['q'];nq=len(q);plane=np.pad(q,((0,0),(0,1)));transfer=d['transfer'];pairs={tuple(p) for p in q}
        assert len(pairs)==len(q) and (0.,0.) not in pairs and all(tuple(-p) not in pairs for p in q)
        for record in r['stages']:
            stage=record['stage'];rng=np.random.default_rng(record['seed']);batch=record['batch'];rotations=Rotation.random(batch,random_state=rng).as_matrix()
            np.testing.assert_array_equal(rotations[0],d[f'stage{stage}_first_rotation'])
            if stage:
                shape=(batch,2,128,nq) if stage==1 else (batch,nq);noise=rng.normal(size=shape)+1j*rng.normal(size=shape)
                np.testing.assert_array_equal(noise[0],d['cal_first_noise' if stage==1 else 'heldout_first_noise'])
            for label,density in [('candidate',rho),('region',rho*mask)]:
                actual=separable(density,plane@rotations[0])*transfer;error=float(np.max(abs(actual-d[f'stage{stage}_{label}_first_direct'])))
                assert error<1e-8;result['maximum_physical_coordinate_difference']=max(result['maximum_physical_coordinate_difference'],error)
        prior_path=BASE/f'candidate-fisher-score-v1/{ds}/summary.json';prior=json.loads(prior_path.read_text());assert prior['complete']
        assert sha(prior_path)==r['input_hashes'][str(prior_path.relative_to(ROOT))]
        prior_array=ROOT/prior['array']['path'];assert sha(prior_array)==prior['array']['sha256']
        with np.load(prior_array) as f:prior_directions={c['key']:f[c['key']+'_direction'] for c in prior['cases']}
        for j,c in enumerate(r['cases']):
            triads=d['triads'] if c['feature_family']=='power_bispectrum' else np.empty((0,3),int);w=d[c['key']+'_direction']
            old=next(x for x in prior['cases'] if x['key']==c['key'])
            assert c['design']==old['design'] and c['baseline_calibration']==old['calibration']
            np.testing.assert_array_equal(w,prior_directions[c['key']]);result['inherited_score_checks']+=1
            threshold=c['design']['threshold']
            # Noise polynomial via four independent feature evaluations and roots.
            mean=d['stage1_candidate_first_direct'];e=d['cal_first_noise'].reshape(-1,nq);amplitudes=np.array([-1.,0.,1.,2.])
            values=np.column_stack([own_score(a*mean+e,triads,w) for a in amplitudes]);coeff=np.linalg.solve(np.vander(amplitudes,4,increasing=True),values.T).T
            events=np.zeros((len(e),16),bool);edges=np.linspace(.9,1.1,17)
            for ell,co in enumerate(coeff):
                p=np.polynomial.Polynomial(co);roots=p.deriv().roots()
                for cell,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
                    points=[lo,hi]+[z.real for z in roots if abs(z.imag)<1e-9 and lo<=z.real<=hi]
                    events[ell,cell]=np.max(p(points))>threshold
            events=events.reshape(2,128,16);np.testing.assert_array_equal(events.mean(axis=1).max(axis=1),d['cal_grouped'][j,0]);np.testing.assert_array_equal(events.max(axis=2).mean(axis=1),d['cal_individual'][j,0]);result['calibration_first_polynomials']+=len(e)
            scores=d['heldout_scores'][j];events=scores>c['design']['threshold'];counts=events.sum(axis=-1);np.testing.assert_array_equal(counts,c['counts']);result['event_count_checks']+=counts.size
            for di,deletion in enumerate(r['deletions']):
                mean=d['stage2_candidate_first_direct']-deletion*d['stage2_region_first_direct']
                for ai,a in enumerate(r['amplitudes']):
                    value=own_score((a*mean+d['heldout_first_noise'])[None],triads,w)[0];error=abs(value-scores[di,ai,0]);assert error<1e-8
                    result['maximum_direct_score_difference']=max(result['maximum_direct_score_difference'],error);result['direct_score_checks']+=1
            # Independently evaluate all seven scalar probability bounds.
            g=d['cal_grouped'][j];raw=d['cal_individual'][j];x=g.mean(axis=1);m=len(x);mu=x.mean()
            assert np.all(g<=raw) and np.all(g>=0) and np.all(raw<=1)
            ordinary=min(1.,mu+radius(x,1.,.001));individual=min(1.,raw.mean()+radius(raw.mean(axis=1),1.,.001))
            rad=radius(x,1.,.00025);jl=max(0.,mu-rad);ju=min(1.,mu+rad)
            product=(g[:,0]-.5)*(g[:,1]-.5);subtract=max(0,jl-.5,.5-ju)**2
            paired=np.clip(product.mean()+radius(product,.5,.0005)-subtract,0,.25)
            sq=(x-.5)**2;unpaired=np.clip(sq.mean()+radius(sq,.25,.0005)-subtract,0,.25)
            rmean=radius(x,1.,.0005);lo=max(0.,mu-rmean);hi=min(1.,mu+rmean);closest=np.clip(.5,lo,hi)
            first=np.sort(x[:m//2]);second=x[m//2:];eps=np.sqrt(np.log(2000)/(2*m));thresholds=np.unique(np.r_[0.,x,1.])
            for bound in c['calibration']['bounds']:
                kappa=bound['kappa'];t=first[max(0,int(np.ceil((1-1/kappa)*len(first)))-1)];excess=np.maximum(0,second-t)
                expected=dict(individual_ratio=min(1.,kappa*individual),grouped_ratio=min(1.,kappa*ordinary),
                    view_variance=min(1.,kappa*ju,ju+np.sqrt((kappa-1)*paired)),
                    mean_only=min(1.,kappa*hi,hi+np.sqrt((kappa-1)*closest*(1-closest))),
                    unpaired_variance=min(1.,kappa*ju,ju+np.sqrt((kappa-1)*unpaired)),
                    cvar_dkw=min(1.,min(t+kappa*(np.maximum(0,x-t).mean()+eps*(1-t)) for t in thresholds)),
                    cvar_split=min(1.,t+kappa*(excess.mean()+radius(excess,1-t,.001))))
                for method,value in expected.items():
                    error=abs(value-bound[method]);assert error<1e-12
                    result['bound_checks']+=1;result['maximum_bound_difference']=max(result['maximum_bound_difference'],error)
            for p in c['projections']+c['repeated_groups']:
                n=p['particles'];k=p['reject_at_count'];bound=p['probability_bound']
                calibration=c['baseline_calibration'] if p['replicas']==32 else c['calibration']
                assert bound==next(b for b in calibration['bounds'] if b['kappa']==p['kappa'])[p['method']]
                assert binom.sf(k-1,n,bound)<=.049+1e-13 and (k==0 or binom.sf(k-2,n,bound)>.049-1e-13);result['critical_value_checks']+=1
                for outcome in p['outcomes']:
                    di=r['deletions'].index(outcome['deletion']);ai=r['amplitudes'].index(outcome['amplitude'])
                    if 'groups' in p:
                        observed=events[di,ai].reshape(-1,n).sum(axis=1);np.testing.assert_array_equal(observed,outcome['event_counts'])
                        assert int(np.sum(observed>=k))==outcome['rejected'];result['group_vector_checks']+=1
                    else:
                        count=int(counts[di,ai]);total=r['heldout_views'];frequency=count/total
                        lo=0. if count==0 else beta.ppf(.025,count,total-count+1);hi=1. if count==total else beta.ppf(.975,count+1,total-count)
                        expected=[binom.sf(k-1,n,frequency),binom.sf(k-1,n,lo),binom.sf(k-1,n,hi)]
                        np.testing.assert_allclose([outcome['rejection_probability'],*outcome['rejection_interval']],expected,rtol=1e-12,atol=1e-14);result['projection_checks']+=1
        print(ds,'inherited scores, allocation bounds and outcome replay complete',flush=True)
    result['complete']=True;args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='input_hashes'}),flush=True)


if __name__=='__main__':main()
