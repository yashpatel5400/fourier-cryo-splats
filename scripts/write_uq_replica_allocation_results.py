#!/usr/bin/env python3
"""Report every allocation outcome without selecting an aggregate test."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development/candidate-replica-allocation-v1'
def writecsv(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def fmt(v):return '<1e-300' if v==0 else f'{v:.4g}'
def main():
    records=[json.loads((BASE/ds/'summary.json').read_text()) for ds in ['10028','10049','10076']];assert all(r['complete'] for r in records)
    v=json.loads((ROOT/'provenance/uncertainty/replica-allocation-independent-verification.json').read_text());assert v['complete']
    rows=[];groups=[];envelopes=[]
    for r in records:
        for c in r['cases']:
            for replicas,cal in [(32,c['baseline_calibration']),(128,c['calibration'])]:
                envelopes.append(dict(dataset=r['dataset'],score=c['key'],replicas=replicas,views=cal['views'],
                    grouped_mean=cal['grouped_mean'],individual_mean=cal['individual_mean'],mean_lower=cal['joint_mean_interval'][0],
                    mean_upper=cal['joint_mean_interval'][1],view_variance_upper=cal['view_variance_upper']))
            for p in c['projections']:
                for o in p['outcomes']:
                    rows.append(dict(dataset=r['dataset'],score=c['key'],replicas=p['replicas'],views=p['calibration_views'],method=p['method'],kappa=p['kappa'],particles=p['particles'],
                        deletion=o['deletion'],amplitude=o['amplitude'],probability_bound=p['probability_bound'],critical_count=p['reject_at_count'],
                        single_image_frequency=o['single_particle_probability'],rejection=o['rejection_probability'],lower=o['rejection_interval'][0],upper=o['rejection_interval'][1]))
            for p in c['repeated_groups']:
                for o in p['outcomes']:
                    groups.append(dict(dataset=r['dataset'],score=c['key'],replicas=p['replicas'],views=p['calibration_views'],method=p['method'],kappa=p['kappa'],particles=p['particles'],
                        deletion=o['deletion'],amplitude=o['amplitude'],critical_count=p['reject_at_count'],groups=p['groups'],rejected=o['rejected'],fraction=o['fraction'],lower=o['pointwise_95_interval'][0],upper=o['pointwise_95_interval'][1]))
    assert len(rows)==37800 and len(groups)==12600
    for name,items in [('projections',rows),('repeated-groups',groups),('envelope-decomposition',envelopes)]:writecsv(BASE/(name+'.csv'),items)
    lines=['# Equal-noise-budget allocation: complete results','',
        'Each allocation uses 2,097,152 conditional noise draws per stack: '
        'M=32,768 views with two groups of 32 versus M=8,192 with two groups '
        'of 128. The larger-group calibration is new; the original one is frozen. '
        'Both are evaluated on the same 131,072 fresh test views, distinct from '
        'the previously viewed Fisher study. Directions and thresholds are unchanged. '
        'This is equal noise-draw count, not equal Fourier calls or wall time, '
        'and is not an optimized allocation or repeated-calibration study.','',
        '## Mean-envelope and outer-uncertainty tradeoff','',
        '| Stack | Score | L | M | Grouped mean | Joint mean upper | Viewing variance upper |',
        '|---|---|---:|---:|---:|---:|---:|']
    for e in envelopes:lines.append(f'| {e["dataset"]} | {e["score"]} | {e["replicas"]} | {e["views"]} | {e["grouped_mean"]:.6f} | {e["mean_upper"]:.6f} | {e["view_variance_upper"]:.6f} |')
    def pick(ds,score,method,kappa,replicas,deletion):return next(x for x in rows if all(x[k]==val for k,val in dict(dataset=ds,score=score,method=method,kappa=kappa,replicas=replicas,deletion=deletion,amplitude=1.,particles=10000).items()))
    lines+=['','## Fixed-score power comparison','',
        'Entries give L=32 / L=128 rejection projections at n=10,000 and '
        'amplitude one. The full CSV retains all seven methods, sample sizes, '
        'caps, amplitudes and correct-null controls. No cross-allocation minimum '
        'is used as a valid test.','',
        '| Stack | Score | Method | kappa | 25% deletion | 50% | Full removal |','|---|---|---|---:|---:|---:|---:|']
    for ds in ['10028','10049','10076']:
        for score in ['power_matched','power_fisher','power_bispectrum_matched','power_bispectrum_fisher']:
            for method in ['view_variance','cvar_split']:
                for kappa in [1.,1.1]:
                    values=[fmt(pick(ds,score,method,kappa,32,d)['rejection'])+' / '+fmt(pick(ds,score,method,kappa,128,d)['rejection']) for d in [.25,.5,1.]]
                    lines.append(f'| {ds} | {score} | {method} | {kappa:g} | '+' | '.join(values)+' |')
    worst=max((r for r in rows if r['deletion']==0 and r['particles']==10000),key=lambda r:r['rejection'])
    g=max((r for r in groups if r['deletion']==0),key=lambda r:r['rejected'])
    lines+=['','## Controls and verification','',
        f'Largest correct-null n=10,000 projection is {worst["rejection"]:.8g} '
        f'[{worst["lower"]:.8g},{worst["upper"]:.8g}] at {worst["dataset"]}, '
        f'{worst["score"]}, L={worst["replicas"]}, {worst["method"]}, '
        f'kappa={worst["kappa"]}, amplitude={worst["amplitude"]}.',
        f'Largest correct-null group count is {g["rejected"]}/128 '
        f'[{g["lower"]:.5g},{g["upper"]:.5g}]. Both calibration realizations are fixed. '
        'All cells share draws and cannot be pooled as independent replications.','',
        f'Independent replay checks {v["inherited_score_checks"]} unchanged score designs, '
        f'{v["bound_checks"]} new bounds, {v["critical_value_checks"]} minimal critical values, '
        f'{v["projection_checks"]} projections and {v["group_vector_checks"]} group vectors. '
        f'It checks {v["calibration_first_polynomials"]} first-view cubics plus direct physical scores. '
        'It does not replay every physical view or provide a certified rounding bound.','',
        'The declared known-simulator assumptions and whole-candidate interpretation '
        'remain unchanged. This comparison does not establish experimental noise/view '
        'calibration, local occupancy coverage, an optimal allocation or a new full-review verdict.','']
    (ROOT/'research/uncertainty/paired-power-v1/REPLICA-ALLOCATION-RESULTS.md').write_text('\n'.join(lines));print('COMPLETE',len(rows),len(groups))
if __name__=='__main__':main()
