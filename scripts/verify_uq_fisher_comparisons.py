#!/usr/bin/env python3
"""Independent conservative difference intervals and joint-event replay."""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from scipy.stats import binom,binomtest
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development/candidate-fisher-score-v1'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=ROOT/'provenance/uncertainty/fisher-comparison-independent-verification.json');args=ap.parse_args()
    if args.output.exists():raise ValueError('Preserve earlier verification')
    path=BASE/'fisher-minus-matched.csv';rows=list(csv.DictReader(path.open()));result=dict(complete=False,input_hashes={str(path.relative_to(ROOT)):sha(path)},joint_event_tables=0,difference_intervals=0,maximum_difference=0.)
    for ds in ['10028','10049','10076']:
        sp=BASE/ds/'summary.json';r=json.loads(sp.read_text());assert r['complete'];p=ROOT/r['array']['path'];assert sha(p)==r['array']['sha256']
        result['input_hashes'][str(sp.relative_to(ROOT))]=sha(sp);result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        with np.load(p) as f:scores=f['heldout_scores']
        cache={}
        for family,j in [('power',0),('power_bispectrum',2)]:
            for di,deletion in enumerate(r['deletions']):
                for ai,amplitude in enumerate(r['amplitudes']):
                    em=scores[j,di,ai]>r['cases'][j]['design']['threshold'];ef=scores[j+1,di,ai]>r['cases'][j+1]['design']['threshold']
                    table=np.bincount(em.astype(int)*2+ef.astype(int),minlength=4);cache[(family,deletion,amplitude)]=table;result['joint_event_tables']+=1
        for row in (x for x in rows if x['dataset']==ds):
            family=row['family'];table=cache[(family,float(row['deletion']),float(row['amplitude']))]
            np.testing.assert_array_equal([int(row[k]) for k in ['neither','fisher_only','matched_only','both']],table)
            total=int(table.sum());cm=int(table[2]+table[3]);cf=int(table[1]+table[3]);n=int(row['particles']);kappa=float(row['kappa'])
            j=0 if family=='power' else 2
            pm,pf=[next(p for p in r['cases'][ix]['projections'] if p['particles']==n and p['kappa']==kappa and p['method']==row['method']) for ix in [j,j+1]]
            km,kf=pm['reject_at_count'],pf['reject_at_count']
            im=binomtest(cm,total).proportion_ci(confidence_level=.975,method='exact');iff=binomtest(cf,total).proportion_ci(confidence_level=.975,method='exact')
            expected=[binom.sf(kf-1,n,cf/total)-binom.sf(km-1,n,cm/total),binom.sf(kf-1,n,iff.low)-binom.sf(km-1,n,im.high),binom.sf(kf-1,n,iff.high)-binom.sf(km-1,n,im.low)]
            error=float(np.max(abs(np.asarray(expected)-[float(row[k]) for k in ['difference','lower','upper']])));assert error<1e-8
            result['maximum_difference']=max(result['maximum_difference'],error);result['difference_intervals']+=1
        print(ds,'joint tables and confidence-difference intervals verified',flush=True)
    assert result['difference_intervals']==9450 and result['joint_event_tables']==90
    result['complete']=True;args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='input_hashes'}),flush=True)
if __name__=='__main__':main()
