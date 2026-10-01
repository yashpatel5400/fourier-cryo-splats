#!/usr/bin/env python3
"""Replay covariance witnesses and independently check selected ridge fields."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from fourier_splats.uq_continuous import continuous_residual_norm
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_paired_covariance import matrix_gaussian_terms
from fourier_splats.uq_covariance_cone_bounds import unrestricted_matrix_growth_upper

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main(output=None):
    record=dict(covariance=[],compression=[],ridge=[],input_hashes={})
    def read(path):
        record['input_hashes'][str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
        return json.loads(path.read_text())
    def arrays(row):
        path=ROOT/row['path'];digest=hashlib.sha256(path.read_bytes()).hexdigest();assert digest==row['sha256']
        record['input_hashes'][row['path']]=digest
        with np.load(path) as f:return {k:f[k].copy() for k in f.files}
    power=read(BASE/'paired-power-enlarged-cone-v1/summary.json')
    compressed=read(BASE/'paired-covariance-finite-view-v1/summary.json')
    full=read(BASE/'paired-covariance-full-frequency-v1/summary.json')
    ridge=read(BASE/'folded-ridge-review3-v1/summary.json')
    for ds in ['10028','10049','10076']:
        p=arrays(power['arrays'][ds]);c=arrays(compressed['arrays'][ds]);f=arrays(full['arrays'][ds])
        transfer=np.sqrt(p['transfer_squared'][0]);ft=p['fourier'][[0,2]]*transfer
        means=np.concatenate([ft.real,ft.imag],axis=-1)
        target=means[0].T@means[0]/len(means[0]);delta=target-means[1].T@means[1]/len(means[1])
        for rank in [8,16,32]:
            basis=c['basis'][:,:rank];ratio=np.linalg.norm(basis.T@delta@basis,'fro')/np.linalg.norm(delta,'fro')
            for index,name in enumerate(['true_map','region_removed']):
                m=c['means'][index,:,:rank];i,j=np.triu_indices(rank)
                features=m[:,i]*m[:,j]*np.where(i==j,1.,np.sqrt(2.))/np.sum(m*m,axis=1)[:,None]
                singular=np.linalg.svd(features,compute_uv=False)
                key=f'{name}_{rank}';weights=c[key+'_coefficients'];approx=np.einsum('n,ni,nj->ij',weights,m,m)
                np.testing.assert_allclose(approx,c[key+'_approximation'],rtol=1e-11,atol=1e-11)
                record['compression'].append(dict(dataset=ds,candidate=name,rank=rank,
                    retained_discrimination_frobenius=float(ratio),symmetric_dimension=rank*(rank+1)//2,
                    feature_rank_relative_tolerance_1e_minus8=int(np.sum(singular>1e-8*singular[0])),
                    smallest_over_largest_singular_value=float(singular[-1]/singular[0])))
        weights=f['coefficients'];assert weights.min()>=0
        approximation=(means[1].T*weights)@means[1]
        np.testing.assert_allclose(approximation,f['approximation'],rtol=1e-10,atol=1e-10)
        fresh=f['fresh_fourier']*transfer;freshmeans=np.concatenate([fresh.real,fresh.imag],axis=1)
        for v in [1.,2.,4.]:
            t=f['union_weights_v'+str(v)];normalizer,u=matrix_gaussian_terms(t)
            growth=float(np.sum(t*target)/v+normalizer)
            normalized=np.sum((freshmeans@u)*freshmeans,axis=1)/np.sum(freshmeans*freshmeans,axis=1)
            assert normalized.max()<1e-9
            old=next(r for r in full['cases'] if r['dataset']==ds and r['candidate']=='region_removed')
            np.testing.assert_allclose(growth,old['growth_after_union_repair'][str(v)]['expected_log_lower'],rtol=1e-10,atol=1e-12)
            upper=unrestricted_matrix_growth_upper(target-approximation,v)['upper']
            record['covariance'].append(dict(dataset=ds,variance_upper=v,union_expected_log=growth,
                maximum_fresh_normalized_constraint=float(normalized.max()),replayed_spectral_upper=upper))
        saved=arrays(ridge['arrays'][ds]);gpath=BASE/f'local-alignment-calibration-v1/{ds}/generators.npz'
        record['input_hashes'][str(gpath.relative_to(ROOT))]=hashlib.sha256(gpath.read_bytes()).hexdigest()
        with np.load(gpath) as generator:
            k,ctf,noise=generator['k'],generator['ctf'],float(generator['noise_std'])
        for name,centers,signs in [('center',[[0,0,0]],[1]),('contrast',[[0,0,.08],[0,0,-.08]],[1,-1])]:
            w=saved[name+'_weights'];exact=continuous_residual_norm(k,ctf,w,noise,centers,signs,.07)
            half=bias_aware_half_width_stable(np.linalg.norm(w),2*exact['residual_norm'],.05)
            old=next(r['fit'] for r in ridge['cases'] if r['dataset']==ds and r['target']==name)
            difference=abs(half-old['half_width'])
            if difference>1e-4:raise ArithmeticError('Direct sinc check differs materially from quadrature width')
            record['ridge'].append(dict(dataset=ds,target=name,independent_sinc_half_width=half,
                stored_quadrature_half_width=old['half_width'],absolute_difference=difference))
        print(ds,'verified',flush=True)
    record['complete']=True
    record['scope']='Numerical replay, SVD diagnostics and independent direct sinc integration, not validated arithmetic.'
    output=ROOT/'provenance/uncertainty/revision5-diagnostic-verification.json' if output is None else Path(output)
    if output.exists():raise RuntimeError('Preserve prior verification')
    output.write_text(json.dumps(record,indent=2)+'\n')
    print('COMPLETE maximum sinc-width difference',max(r['absolute_difference'] for r in record['ridge']),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output');args=parser.parse_args();main(args.output)
