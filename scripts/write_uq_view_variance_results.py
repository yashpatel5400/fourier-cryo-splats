#!/usr/bin/env python3
"""All-case report for a selected, explicitly post-outcome follow-up."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    records=[json.loads((BASE/f'bispectrum-view-variance-v1/{ds}/summary.json').read_text()) for ds in ['10028','10049']]
    assert all(r['complete'] for r in records)
    verified=json.loads((ROOT/'provenance/uncertainty/view-variance-independent-verification.json').read_text());assert verified['complete']
    lines=['# Conditional-noise and viewing-variance refinement: complete results','',
        'This is a selected post-outcome follow-up of all four nonzero ranged-amplitude '
        'directions on 10028/10049. It is not a new three-stack positive result; the '
        'fixed-amplitude and zero 10076 directions remain in the previous inventory. '
        'Scores and thresholds are unchanged. Each stack has 32,768 new Haar '
        'calibration views, two groups of 32 independent noises per view, and '
        '131,072 new independent Haar held-out views. Both true and removed candidates '
        'and all three matched-simulation-budget procedures are retained.','',
        'The three methods use the individual-noise amplitude envelope, grouped '
        'amplitude envelope, and grouped viewing-variance bound, respectively. '
        'All use the same assumed [.9,1.1] amplitude interval, known CTF/proper '
        'Gaussian noise, and supplied viewing-density ratio. The error argument '
        'uses independent views as its sample units. Conditional replicas are not '
        'counted as independent views. No experimental parameter is calibrated here.','',
        '## Removed-candidate rejection projections at n=10,000','',
        'Each entry is a binomial projection on the Haar true-map alternative '
        'with a pointwise 95% interval from held-out event frequencies. These are '
        'not observed repeated-dataset power or simultaneous confidence intervals. '
        'Changing kappa changes the null guarantee and threshold, not the '
        'held-out viewing distribution.','',
        '| Stack | Score | Procedure | kappa=1 | kappa=1.1 | kappa=2 | kappa=5 |',
        '|---|---|---|---:|---:|---:|---:|']
    rows=[]
    for r in records:
        ds=r['dataset']
        for c in r['cases']:
            for method in ['individual_ratio','grouped_ratio','view_variance']:
                displays=[]
                for kappa in [1.,1.1,2.,5.]:
                    p=next(p for p in c['projections'] if p['candidate']=='removed' and p['method']==method and p['kappa']==kappa and p['particles']==10000)
                    x=p['true'];lo,hi=x['rejection_interval'];displays.append(f"{x['rejection_probability']:.3g} [{lo:.3g}, {hi:.3g}]")
                lines.append(f"| {ds} | {c['key']} | {method} | "+' | '.join(displays)+' |')
            for p in c['projections']:
                row=dict(dataset=ds,score=c['key'],candidate=p['candidate'],procedure=p['method'],
                    kappa=p['kappa'],particles=p['particles'],probability_bound=p['probability_bound'],
                    reject_at_count=p['reject_at_count'])
                for label in ['true','removed']:
                    row[label+'_frequency']=p[label]['single_particle_probability']
                    row[label+'_rejection']=p[label]['rejection_probability']
                    row[label+'_lower'],row[label+'_upper']=p[label]['rejection_interval']
                rows.append(row)
    maximum=max(c['maximum_coordinate_difference'] for r in records for s in r['stages'] for c in s['operator_checks'])
    lines+=['',f'The full CSV retains all {len(rows)} projections at n=1,000/10,000/100,000 '
        'and kappa=1/1.01/1.1/2/5, including both correct-null controls. '
        f"The two runs take {records[0]['seconds']:.2f} and {records[1]['seconds']:.2f} seconds "
        'in overlapping Mac processes; their sum is not elapsed wall time. '
        f'Eight first-view physical-cell sums agree within {maximum:.3g} per coordinate. '
        f"Independent replay checks {verified['bound_checks']} probability bounds, "
        f"{verified['critical_value_checks']} critical values, all eight held-out count arrays "
        'and 512 first-view noise polynomials across 16 amplitude cells.','',
        '## Interpretation and limitations','',
        'The grouping reduces noise-dependent amplitude selection in the simulation '
        'envelope. The paired-group product estimates variation of its conditional '
        'mean across views, enabling a less crude density-ratio penalty. Its '
        'ingredients are classical empirical Bernstein concentration, conditional '
        'covariance and chi-squared/Cauchy--Schwarz robustness, with primary '
        'reading scopes documented separately. A numerical improvement is not '
        'a novelty verdict. The score still knows an oracle alternative, only '
        'one transfer profile and full-removal alternative are tested, and all '
        'viewing laws here are simulated Haar. Experimental noise/view/amplitude '
        'calibration, heterogeneity, smaller changes and useful learned candidates '
        'remain unresolved. All three full Fable reviews remain rejections.','']
    output=BASE/'bispectrum-view-variance-v1/projections.csv'
    with output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    (ROOT/'research/uncertainty/paired-power-v1/VIEW-VARIANCE-RESULTS.md').write_text('\n'.join(lines))
    print('Wrote',len(rows),'complete projections')


if __name__=='__main__':main()
