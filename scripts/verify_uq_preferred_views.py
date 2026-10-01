#!/usr/bin/env python3
"""Replay scores, event counts, thresholds and sampling metadata independently."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import binom
from fourier_splats.uq_bispectrum import moment_features
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'provenance/uncertainty/preferred-views-independent-verification.json');args=parser.parse_args()
    out=args.output
    if out.exists():raise ValueError('Preserve previous verification')
    result=dict(complete=False,input_hashes={},geometry=[],first_score_checks=0,event_count_checks=0,critical_value_checks=0,
        repeated_group_checks=0,maximum_first_score_error=0.)
    h=json.loads((BASE/'bispectrum-hull-v1/summary.json').read_text());original=json.loads((BASE/'paired-power-enlarged-cone-v1/summary.json').read_text())
    for ds in ['10028','10049']:
        sp=BASE/f'bispectrum-preferred-view-v1/{ds}/summary.json';r=json.loads(sp.read_text());assert r['complete']
        cp=BASE/f'bispectrum-view-variance-v1/{ds}/summary.json';cal=json.loads(cp.read_text());assert cal['complete']
        dp=ROOT/h['arrays'][ds]['path'];assert sha(dp)==h['arrays'][ds]['sha256']
        op=ROOT/original['arrays'][ds]['path'];assert sha(op)==original['arrays'][ds]['sha256']
        for p in [sp,cp,dp,op]:result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        with np.load(dp) as f:d={key:f[key].copy() for key in f.files}
        q=d['q'];assert q.shape==(220,2) and len(np.unique(q,axis=0))==len(q) and np.all(np.linalg.norm(q,axis=1)>0)
        pairs={tuple(x) for x in q};assert all(tuple(-x) not in pairs for x in q)
        assert d['transfer'].shape==(len(q),)
        with np.load(op) as f:np.testing.assert_array_equal(d['transfer'],np.sqrt(f['transfer_squared'][0]))
        result['geometry'].append(dict(dataset=ds,frequencies=len(q),nonzero_unique_no_antipodes=True,first_transfer_profile_exactly_matches=True))
        for law in r['laws']:
            ap=ROOT/law['array']['path'];assert sha(ap)==law['array']['sha256'];result['input_hashes'][str(ap.relative_to(ROOT))]=sha(ap)
            with np.load(ap) as f:data={key:f[key].copy() for key in f.files}
            # Reproduce the first rejection-sampling batch without its helper.
            rng=np.random.default_rng(law['seed']);kept=[];remaining=r['batch']
            while remaining:
                raw=Rotation.random(max(16,int(np.ceil(remaining*law['kappa']*1.05))),random_state=rng).as_matrix()
                accepted=raw[np.abs(raw[:,2,law['axis']])>=law['cut']][:remaining]
                kept.append(accepted);remaining-=len(accepted)
            rotations=np.concatenate(kept);noise=rng.normal(size=(r['batch'],len(q)))+1j*rng.normal(size=(r['batch'],len(q)))
            np.testing.assert_array_equal(rotations[0],data['example_rotation']);np.testing.assert_array_equal(noise[0],data['example_noise'])
            assert law['minimum_absolute_coordinate']>=law['cut'] and law['used']==r['total_per_law']
            for c in law['cases']:
                j=next(j for j,z in enumerate(cal['cases']) if z['key']==c['key']);a_index=r['amplitudes'].index(c['amplitude'])
                frozen=cal['cases'][j];assert c['threshold']==frozen['threshold']
                triads=d['triads'] if c['key'].startswith('power_bispectrum') else np.empty((0,3),int)
                for label in ['true','removed']:
                    mean=data[label+'_example_direct'];w=d[c['key']+'_direction'];a=c['amplitude']
                    score=float(moment_features(a*mean+noise[0],triads,1.)@w)
                    error=abs(score-data[label][j,a_index,0]);assert error<1e-8
                    result['maximum_first_score_error']=max(result['maximum_first_score_error'],error);result['first_score_checks']+=1
                    events=data[label][j,a_index]>c['threshold'];assert int(events.sum())==c['counts'][label];result['event_count_checks']+=1
                    for p in c['repeated_groups']:
                        counts=events.reshape(-1,p['particles']).sum(axis=1)
                        np.testing.assert_array_equal(counts,p[label]['event_counts'])
                        assert int(np.sum(counts>=p['reject_at_count']))==p[label]['rejected'];result['repeated_group_checks']+=1
                for p in c['projections']+c['repeated_groups']:
                    archived=next(b for b in frozen['calibration'][p['candidate']]['bounds'] if b['kappa']==law['kappa'])[p['method']]
                    assert archived==p['probability_bound'];n=p['particles'];k=p['reject_at_count']
                    assert binom.sf(k-1,n,archived)<=.049+1e-13
                    assert k==0 or binom.sf(k-2,n,archived)>.049-1e-13
                    result['critical_value_checks']+=1
        print(ds,'verified',flush=True)
    result['complete']=True;out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='input_hashes'}),flush=True)


if __name__=='__main__':main()
