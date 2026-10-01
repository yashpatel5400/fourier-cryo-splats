#!/usr/bin/env python3
"""Complete preferred-view and classical-comparator reports, with all controls."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def csv_write(path,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)


def fmt(x):return '<1e-300' if x==0 else f'{x:.4g}'


def main():
    verified=json.loads((ROOT/'provenance/uncertainty/preferred-views-independent-verification.json').read_text());assert verified['complete']
    records=[json.loads((BASE/f'bispectrum-preferred-view-v1/{ds}/summary.json').read_text()) for ds in ['10028','10049']]
    assert all(r['complete'] for r in records)
    rows=[];replicates=[];table=[];diagnostics=[]
    for r in records:
        for law in r['laws']:
            diagnostics.append(dict(dataset=r['dataset'],axis=law['axis'],kappa=law['kappa'],
                observed_acceptance=law['proposal_acceptance_fraction'],expected_acceptance=1/law['kappa'],
                second_moment=law['coordinate_square_mean'],expected_second_moment=law['theoretical_coordinate_square_mean']))
            for c in law['cases']:
                identity=dict(dataset=r['dataset'],axis=law['axis'],kappa=law['kappa'],score=c['key'],amplitude=c['amplitude'])
                for p in c['projections']:
                    row=dict(**identity,candidate=p['candidate'],method=p['method'],particles=p['particles'],probability_bound=p['probability_bound'],reject_at_count=p['reject_at_count'])
                    for label in ['true','removed']:
                        row[label+'_frequency']=p[label]['single_particle_probability'];row[label+'_rejection']=p[label]['rejection_probability']
                        row[label+'_lower'],row[label+'_upper']=p[label]['rejection_interval']
                    rows.append(row)
                for p in c['repeated_groups']:
                    row=dict(**identity,candidate=p['candidate'],method=p['method'],particles=p['particles'],groups=p['groups'],reject_at_count=p['reject_at_count'])
                    for label in ['true','removed']:
                        row[label+'_rejected']=p[label]['rejected'];row[label+'_fraction']=p[label]['fraction']
                        row[label+'_lower'],row[label+'_upper']=p[label]['pointwise_95_interval']
                    replicates.append(row)
    lines=['# Fresh preferred-view controls: complete results','',
        'All 18 law/stack runs finish (nine laws on each of two stacks). Each law has '
        '65,536 new independent views, with density kappa on an axis-specific '
        'Haar cap union. There are 1,179,648 independently drawn rotations in total; '
        'maps, amplitudes and scores share draws. The physical amplitude is independent '
        'of the measurement noise in this generator. All cases in the frozen protocol '
        'are retained, including both correct-null maps and all .9/1/1.1 amplitudes.','',
        '## Removed-candidate rejection on true-map data at n=10,000','',
        'The table reports the range of point projections across all three spatial '
        'axes, using the paired view-variance procedure. These ranges are not '
        'confidence intervals or power guarantees over all viewing laws. Pointwise '
        '95% intervals and every method are in the complete CSV. Actual rejected '
        'groups of size 1,024 are reported separately.','',
        '| Stack | Score | kappa | a=.9 | a=1 | a=1.1 |','|---|---|---:|---:|---:|---:|']
    for ds in ['10028','10049']:
        for score in ['power_range09_11','power_bispectrum_range09_11']:
            for kappa in [1.1,2.,5.]:
                displays=[]
                for a in [.9,1.,1.1]:
                    subset=[r for r in rows if r['dataset']==ds and r['score']==score and r['kappa']==kappa and r['amplitude']==a and r['candidate']=='removed' and r['method']=='view_variance' and r['particles']==10000]
                    values=[r['true_rejection'] for r in subset];displays.append(f'{fmt(min(values))}--{fmt(max(values))}')
                lines.append(f'| {ds} | {score} | {kappa:g} | '+' | '.join(displays)+' |')
    nulls=[r[r['candidate']+'_rejection'] for r in rows if r['particles']==10000]
    null_groups=[r[r['candidate']+'_rejected'] for r in replicates]
    maximum=max(r['operator_checks'][j]['maximum_coordinate_difference'] for s in records for r in s['laws'] for j in range(2))
    lines+=['','## Correct-null and numerical controls','',
        f'Across all {len(replicates)} correct-null method/law/amplitude/score cells, '
        f'the largest observed rejection count is {max(null_groups)}/64 groups of 1,024 '
        'particles. The associated exact 95% interval is wide; absence of an observed '
        'excess is not a proof of unconditional 5% control. Calibration was held fixed, '
        'and cells share draws, so pooling these outcomes as independent trials would '
        'be wrong. At n=10,000, the largest correct-null point projection across all '
        f'procedures is {max(nulls):.6g}.','',
        f'The full CSV contains {len(rows)} binomial projections and a separate CSV '
        f'contains {len(replicates)} pairs of actual group-rejection counts. '
        'Every 64-vector of observed counts remains in JSON. '
        f'The independent verifier replays {verified["event_count_checks"]} event counts, '
        f'{verified["critical_value_checks"]} critical values and {verified["repeated_group_checks"]} '
        'group-count vectors; the first 216 direct scores agree within '
        f'{verified["maximum_first_score_error"]:.3g}. Thirty-six physical-cell first-view '
        f'checks agree within {maximum:.3g} per Fourier coordinate. Frequencies are '
        'unique, nonzero and have no antipodal pairs, and the first transfer profile '
        'matches its saved source exactly. Only these declared score/operator checks '
        'were independently replayed; this is not interval arithmetic.','',
        f'The two overlapping Mac runs take {records[0]["seconds"]:.2f} and '
        f'{records[1]["seconds"]:.2f} seconds. Their sum is not elapsed wall time.','',
        '## Interpretation','',
        'The prescribed preferential viewing does not trigger an observed type-I '
        'failure in these controls. Power changes materially with amplitude, so '
        'amplitude-one Haar results do not summarize this experiment. The tests remain '
        'weak under strong view concentration. Axis caps are not the extremal '
        'event-probability tilt; an adversarial-law study remains outstanding. The '
        'known proper noise, exact transfer profile, supplied viewing-density cap '
        'and oracle full-removal score are unchanged limitations. No experimental '
        'viewing law, uncertainty coverage or three-stack success is established.','']
    prefix=BASE/'bispectrum-preferred-view-v1';csv_write(prefix/'projections.csv',rows);csv_write(prefix/'repeated-groups.csv',replicates);csv_write(prefix/'sampler-diagnostics.csv',diagnostics)
    (ROOT/'research/uncertainty/paired-power-v1/PREFERRED-VIEW-RESULTS.md').write_text('\n'.join(lines))
    risk=json.loads((BASE/'bispectrum-view-risk-v1/summary.json').read_text());assert risk['complete']
    lines=['# Classical risk comparators and Monte Carlo slack: complete results','',
        'Post-outcome analysis of the same two calibration archives, motivated by '
        'focused audit V1--V5. All eight score/candidate cases are retained. The '
        'grouped methods require amplitude conditionally independent of noise given '
        'view. The individual envelope allows a larger null. None of the four new '
        'procedures is selected by taking an unadjusted minimum with another.','',
        '## Removed-candidate probability bounds at kappa=1.1','',
        '| Stack | Score | Paired variance | CVaR DKW | CVaR split EB | Mean only | Unpaired variance |',
        '|---|---|---:|---:|---:|---:|---:|']
    riskrows=[];slack=[]
    for c in risk['cases']:
        ds=c['dataset'];old=json.loads((BASE/f'bispectrum-view-variance-v1/{ds}/summary.json').read_text())
        original=next(z for z in old['cases'] if z['key']==c['key'])
        if c['candidate']=='removed':
            b=next(z for z in c['baselines']['bounds'] if z['kappa']==1.1)
            paired=next(z for z in original['calibration']['removed']['bounds'] if z['kappa']==1.1)['view_variance']
            lines.append(f'| {ds} | {c["key"]} | {paired:.6f} | '+' | '.join(f'{b[k]:.6f}' for k in ['cvar_dkw','cvar_split','mean_only','unpaired_variance'])+' |')
        slack.append(dict(dataset=ds,score=c['key'],candidate=c['candidate'],**c['decomposition']))
        for p in c['projections']:
            row=dict(dataset=ds,score=c['key'],candidate=c['candidate'],method=p['method'],kappa=p['kappa'],particles=p['particles'],
                probability_bound=p['probability_bound'],reject_at_count=p['reject_at_count'],critical_minus_paired=p['critical_minus_paired'])
            for label in ['true','removed']:
                row[label+'_rejection']=p[label]['rejection_probability'];row[label+'_lower'],row[label+'_upper']=p[label]['rejection_interval']
            riskrows.append(row)
    lines+=['','## Variance-bound decomposition (all candidates)','',
        '| Stack | Score | Candidate | Raw product mean | Square-root term | Range term | Distance subtraction | Final upper V |',
        '|---|---|---|---:|---:|---:|---:|---:|']
    for s in slack:
        lines.append('| '+f'{s["dataset"]} | {s["score"]} | {s["candidate"]} | '+' | '.join(f'{s[k]:.6g}' for k in ['product_estimate','square_root_term','range_term','distance_subtraction','variance_upper'])+' |')
    lines+=['',f'All {len(riskrows)} new projections and paired critical-count differences '
        'are in CSV. For removed 10049 candidates, raw centered-product estimates '
        'are negative; adding the finite Monte Carlo concentration terms makes the '
        'upper bounds positive. These results confirm slack-dominated variance '
        'bounds and do not measure negligible population view variance.','',
        'The paired method gives smaller bounds than these implemented CVaR and '
        'unpaired comparators at kappa=1.1, for this particular budget. This does '
        'not refute exact population CVaR optimality: CVaR is applied to noisy '
        'view-group means, whereas the paired product targets their conditional '
        'mean variance. Nor does it prove superiority under optimized allocation, '
        'other replication budgets, or independent calibration repetitions. '
        'Those comparisons remain unperformed. At fixed saved frequencies, '
        'method ordering follows the critical counts directly; overlapping '
        'single-image probability intervals do not test a paired method difference.','',
        'Six targeted tests check the CVaR density linear program (including '
        'fractional boundary atoms), explicit DKW breakpoint evaluation, conditional '
        'Jensen ordering, a noise-adaptive-amplitude counterexample, near-degenerate '
        'cubics and the variance decomposition. They do not prove scientific novelty. '
        'The primary reading scope and all formulas are in the frozen protocol.','']
    csv_write(BASE/'bispectrum-view-risk-v1/projections.csv',riskrows);csv_write(BASE/'bispectrum-view-risk-v1/variance-decomposition.csv',slack)
    (ROOT/'research/uncertainty/paired-power-v1/VIEW-RISK-BASELINES-RESULTS.md').write_text('\n'.join(lines))
    print('Wrote',len(rows),'preferred projections,',len(replicates),'replicate cells,',len(riskrows),'risk projections')


if __name__=='__main__':main()
