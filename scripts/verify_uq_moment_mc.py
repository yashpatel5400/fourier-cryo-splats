#!/usr/bin/env python3
"""Independently replay event counts, first scores and binomial inequalities."""
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
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,
        default=ROOT/'provenance/uncertainty/moment-mc-independent-verification.json');args=parser.parse_args()
    if args.output.exists():raise ValueError('Preserve previous verification')
    sp=BASE/'bispectrum-monte-carlo-v1/summary.json';r=json.loads(sp.read_text());assert r['complete']
    hp=BASE/'bispectrum-hull-v1/summary.json';h=json.loads(hp.read_text())
    result=dict(complete=False,input_hashes={str(p.relative_to(ROOT)):sha(p) for p in [sp,hp]},
        event_count_checks=0,critical_value_checks=0,first_score_checks=[],calibration_tail_residuals=[],
        counts_changed_by_1e8_threshold_perturbation=0,minimum_threshold_distance=None)
    distances=[]
    for ds,record in r['arrays'].items():
        p=ROOT/record['path'];assert sha(p)==record['sha256']
        ap=ROOT/h['arrays'][ds]['path'];assert sha(ap)==h['arrays'][ds]['sha256']
        for path in [p,ap]:result['input_hashes'][str(path.relative_to(ROOT))]=sha(path)
        with np.load(p) as f:archive={k:f[k].copy() for k in f.files}
        with np.load(ap) as f:data={k:f[k].copy() for k in f.files}
        cases=[c for c in r['cases'] if c['dataset']==ds and 'skipped' not in c]
        for stage in [0,1]:
            stage_record=next(s for s in r['stages'] if s['dataset']==ds and s['stage']==stage)
            rng=np.random.default_rng(stage_record['seed']);rotations=Rotation.random(r['batch'],random_state=rng).as_matrix()
            np.testing.assert_allclose(rotations[0],archive[f'stage{stage}_example_rotation'],atol=0,rtol=0)
            n=len(data['q']);noise=(rng.normal(size=(r['batch'],n))+1j*rng.normal(size=(r['batch'],n)))[0]
            for j,c in enumerate(cases):
                key=c['key'];w=data[key+'_direction']
                triads=data['triads'] if key.startswith('power_bispectrum') else np.empty((0,3),int)
                for label in ['true','removed']:
                    mean=archive[f'stage{stage}_{label}_example_direct']
                    # Independent polynomial interpolation and companion roots,
                    # not the production coefficient/product/root calculation.
                    amplitudes=np.array([-1.,0.,1.,2.])
                    values=moment_features(amplitudes[:,None]*mean+noise,triads,1.)@w
                    coefficients=np.linalg.solve(np.vander(amplitudes,4,increasing=True),values)
                    polynomial=np.polynomial.Polynomial(coefficients)
                    lo,hi=c['amplitude_interval'];points=[lo,hi]
                    for root in polynomial.deriv().roots():
                        if abs(root.imag)<1e-9 and lo<=root.real<=hi:points.append(float(root.real))
                    actual=float(max(polynomial(points)));saved=archive[f'stage{stage}_{label}_envelope'][j,0]
                    fixed=float(moment_features(mean+noise,triads,1.)@w)
                    np.testing.assert_allclose([actual,fixed],
                        [saved,archive[f'stage{stage}_{label}_fixed'][j,0]],atol=2e-8,rtol=1e-8)
                    result['first_score_checks'].append(dict(dataset=ds,stage=stage,key=key,map=label,
                        envelope_discrepancy=float(abs(actual-saved))))
        for j,c in enumerate(cases):
            threshold=c['threshold']
            for name,count in c['counts'].items():
                scores=archive[name][j];actual=int(np.sum(scores>threshold));assert actual==count['count']
                for adjustment in [-1e-8,1e-8]:
                    result['counts_changed_by_1e8_threshold_perturbation']+=abs(int(np.sum(scores>threshold+adjustment))-actual)
                distances.append(float(np.min(abs(scores-threshold))));result['event_count_checks']+=1
            for candidate,cal in c['calibration'].items():
                residual=abs(float(binom.cdf(cal['event_count'],r['samples_per_stage'],cal['probability_upper']))-.001)
                assert residual<2e-10;result['calibration_tail_residuals'].append(residual)
            for p in c['projections']:
                n=p['particles'];k=p['reject_at_count'];prob=p['calibrated_probability_bound']
                assert binom.sf(k-1,n,prob)<=.049+1e-13
                assert k==0 or binom.sf(k-2,n,prob)>.049-1e-13
                result['critical_value_checks']+=1
        print(ds,'all event counts, independent first-score polynomials and critical values replayed',flush=True)
    result.update(complete=True,minimum_threshold_distance=min(distances))
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('COMPLETE',result['event_count_checks'],result['critical_value_checks'],
        len(result['first_score_checks']),result['counts_changed_by_1e8_threshold_perturbation'],flush=True)


if __name__=='__main__':main()
