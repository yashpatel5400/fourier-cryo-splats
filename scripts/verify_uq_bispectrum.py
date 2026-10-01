#!/usr/bin/env python3
"""Replay mixtures and independently check contrast variance by complex moments."""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
from fourier_splats.uq_bispectrum import moment_features,MomentContrast

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


class IndependentComplexMoments:
    """Enumerate polynomial monomials; independent of real derivative code."""
    def __init__(self,n,triads,w):
        t=len(triads);terms=[];weights=[]
        for j,x in enumerate(w[:n]/2):
            if x:terms.append({j:(1,1)});weights.append(complex(x))
        b=(w[n:n+t]-1j*w[n+t:])/4
        for (a,k,c),x in zip(triads,b):
            if not x:continue
            terms.append({a:(1,0),k:(1,0),c:(0,1)});weights.append(x)
            terms.append({a:(0,1),k:(0,1),c:(1,0)});weights.append(x.conjugate())
        def encoded(items):
            indices=np.zeros((len(items),6),int);pp=indices.copy();qq=indices.copy()
            for row,term in enumerate(items):
                for column,(j,(p,q)) in enumerate(term.items()):
                    indices[row,column]=j;pp[row,column]=p;qq[row,column]=q
            return indices,pp,qq
        self.single=encoded(terms);pairs=[];s_indices=[];t_indices=[];factors=[]
        for s,first in enumerate(terms):
            for t in range(s,len(terms)):
                second=terms[t]
                if not (first.keys()&second.keys()):continue
                merged={}
                for j in first.keys()|second.keys():
                    p,q=first.get(j,(0,0));u,v=second.get(j,(0,0));merged[j]=(p+u,q+v)
                pairs.append(merged);s_indices.append(s);t_indices.append(t)
                factors.append(weights[s]*weights[t]*(1 if s==t else 2))
        self.pairs=encoded(pairs);self.left=np.array(s_indices);self.right=np.array(t_indices)
        self.factors=np.array(factors)

    @staticmethod
    def evaluate(encoded,mean,v):
        index,p,q=encoded;m=mean[index];values=np.zeros_like(m)
        # E[Y^p conj(Y)^q] for proper complex Gaussian with E|noise|²=2v.
        for a in range(3):
            for b in range(3):
                selected=(p==a)&(q==b)
                if not np.any(selected):continue
                z=m[selected];value=np.zeros(len(z),complex)
                for k in range(min(a,b)+1):
                    value+=math.comb(a,k)*math.comb(b,k)*math.factorial(k)*(2*v)**k*z**(a-k)*z.conj()**(b-k)
                values[selected]=value
        return np.prod(values,axis=1)

    def variance(self,mean,v):
        single=self.evaluate(self.single,mean,v)
        joint=self.evaluate(self.pairs,mean,v)
        total=np.sum(self.factors*(joint-single[self.left]*single[self.right]))
        if abs(total.imag)>1e-8*max(1.,abs(total.real)):
            raise ArithmeticError('A real contrast has complex variance')
        return float(total.real)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,
        default=ROOT/'provenance/uncertainty/bispectrum-independent-verification.json');args=parser.parse_args()
    if args.output.exists():raise ValueError('Choose an unused output path; preserve prior verification records')
    result=dict(complete=False,cases=[],input_hashes={},scope='Saved-array numerical replay and independent complex Gaussian moment enumeration, not calibrated arithmetic or continuous-orbit proof.')
    sp=BASE/'bispectrum-hull-v1/summary.json';r=json.loads(sp.read_text())
    pp=BASE/'bispectrum-adversarial-pose-v1/summary.json';p=json.loads(pp.read_text())
    op=BASE/'paired-power-enlarged-cone-v1/summary.json';old=json.loads(op.read_text())
    for path in [sp,pp,op]:result['input_hashes'][str(path.relative_to(ROOT))]=sha(path)
    assert r['complete'] and p['complete']
    for ds in ['10028','10049','10076']:
        files=[]
        for record in [r,p,old]:
            path=ROOT/record['arrays'][ds]['path'];assert sha(path)==record['arrays'][ds]['sha256'];files.append(path)
            result['input_hashes'][str(path.relative_to(ROOT))]=sha(path)
        data,poses,original=[dict(np.load(path)) for path in files]
        means=original['fourier'][[0,2]]*data['transfer']
        for case in r['cases']:
            if case['dataset']!=ds or case['candidate']!='region_removed':continue
            key=case['features']+('_fixed' if len(case['amplitudes'])==1 else '_range09_11')
            triads=data['triads'] if case['features']=='power_bispectrum' else np.empty((0,3),int)
            target=moment_features(means[0],triads).mean(axis=0)
            atoms=np.concatenate([moment_features(a*means[1],triads) for a in case['amplitudes']])
            w=data[key+'_weights'];assert w.min()>=0 and abs(w.sum()-1)<1e-12
            approximation=w@atoms;direction=data[key+'_direction']
            upper=float(np.linalg.norm(target-approximation))
            lower=max(0.,float(direction@target-np.max(atoms@direction)))
            np.testing.assert_allclose(upper,case['distance_upper'],rtol=1e-10,atol=1e-11)
            np.testing.assert_allclose(lower,case['distance_lower'],rtol=1e-10,atol=1e-11)
            row=dict(dataset=ds,key=key,distance_upper=upper,distance_lower=lower,variance_checks=[])
            if np.linalg.norm(direction):
                independent=IndependentComplexMoments(len(data['q']),triads,direction)
                contrast=MomentContrast(len(data['q']),triads,direction)
                worst=next(c['worst_start'] for c in p['cases'] if c['dataset']==ds and c['key']==key)
                for label,mean in [('true_0',means[0,0]),('true_143',means[0,143]),('adversarial',poses[key+'_means'][worst])]:
                    for amplitude in [.9,1.1]:
                        m=amplitude*mean;derivative=float(contrast.mean_variance(m)[1][0]);wick=independent.variance(m,1.)
                        np.testing.assert_allclose(derivative,wick,rtol=2e-10,atol=2e-10)
                        row['variance_checks'].append(dict(mean=label,amplitude=amplitude,
                            derivative_variance=derivative,complex_moment_variance=wick,absolute_difference=abs(derivative-wick)))
            result['cases'].append(row);print(ds,key,'replayed; variance checks',len(row['variance_checks']),flush=True)
    result['complete']=True;result['variance_comparisons']=sum(len(c['variance_checks']) for c in result['cases'])
    result['maximum_variance_difference']=max(x['absolute_difference'] for c in result['cases'] for x in c['variance_checks'])
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('COMPLETE',len(result['cases']),result['variance_comparisons'],result['maximum_variance_difference'],flush=True)


if __name__=='__main__':main()
