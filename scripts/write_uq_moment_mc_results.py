#!/usr/bin/env python3
"""Write every Monte Carlo projection and a clearly limited compact report."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    r=json.loads((BASE/'bispectrum-monte-carlo-v1/summary.json').read_text());assert r['complete']
    verification=json.loads((ROOT/'provenance/uncertainty/moment-mc-independent-verification.json').read_text())
    assert verification['complete']
    rows=[];lines=['# Monte Carlo candidate validation: bounded-view feasibility results','',
        f"The frozen study completes in {r['seconds']:.3f} seconds. Each nonzero stack uses "
        '131,072 independent continuous Haar calibration draws and 131,072 held-out draws. '
        'The true and region-removed maps share rotations/noise within a stage. '
        'Calibration and held-out stages are independent. Four zero 10076 directions '
        'are retained as uninformative, with no projections or random draws. '
        'The study has no new experimental images or reconstruction results.','',
        'The test assumes a known proper Gaussian noise/CTF model, an amplitude interval '
        'and a specified bound kappa on the viewing density relative to Haar measure. '
        'Its alpha=.05 guarantee includes delta=.001 calibration failure; it does '
        'not condition on arbitrary realized poses or establish those assumptions '
        'from experimental data. Each candidate/test is assessed separately.','',
        '## All nonzero contrasts at 100,000 tested particles','',
        'Entries are projected rejection probabilities on the true-map Haar alternative, '
        'using the independently calibrated region-removed candidate. Brackets are '
        'pointwise 95% intervals obtained from held-out single-particle probabilities. '
        'They are not observed repeated-dataset power or simultaneous intervals.','',
        '| Stack | Contrast | kappa=1 | kappa=1.01 | kappa=1.1 | kappa=2 |',
        '|---|---|---:|---:|---:|---:|']
    for c in r['cases']:
        if 'skipped' in c:continue
        display=[]
        for p in c['projections']:
            row=dict(dataset=c['dataset'],contrast=c['key'],candidate=p['candidate'],
                kappa=p['kappa'],particles=p['particles'],threshold=c['threshold'],
                calibrated_probability_bound=p['calibrated_probability_bound'],reject_at_count=p['reject_at_count'])
            for label in ['true_fixed','removed_fixed',p['candidate']+'_envelope']:
                x=p[label];row[label+'_event_frequency']=x['single_particle_probability']
                row[label+'_rejection_probability']=x['rejection_probability']
                row[label+'_rejection_lower']=x['rejection_interval'][0]
                row[label+'_rejection_upper']=x['rejection_interval'][1]
            rows.append(row)
            if p['candidate']=='removed' and p['particles']==100000:
                x=p['true_fixed'];lo,hi=x['rejection_interval']
                display.append(f"{x['rejection_probability']:.3g} [{lo:.3g}, {hi:.3g}]")
        lines.append(f"| {c['dataset']} | {c['key']} | "+' | '.join(display)+' |')
    out=BASE/'bispectrum-monte-carlo-v1/projections.csv'
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with out.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys,lineterminator='\n');writer.writeheader();writer.writerows(rows)
    true_controls=[p['true_fixed']['rejection_probability'] for c in r['cases'] if 'projections' in c
        for p in c['projections'] if p['candidate']=='true']
    removed_controls=[p['removed_fixed']['rejection_probability'] for c in r['cases'] if 'projections' in c
        for p in c['projections'] if p['candidate']=='removed']
    maximum=max(x['maximum_coordinate_difference'] for s in r['stages'] for x in s['operator_checks'])
    lines.extend(['',f'The full CSV retains all {len(rows)} projections at n=1,000, 10,000 and 100,000, '
        'including both correct-map and removed-map null controls. Maximum point '
        f'projections of null rejection are {max(true_controls):.6g} for the true-map '
        f'control and {max(removed_controls):.6g} for removed-map amplitude-one data; '
        'their complete intervals remain in the CSV. These are projections conditional '
        'on the realized calibration, not empirical unconditional false-positive rates.','',
        f'All eight independent first-view physical-cell sums agree with the NUFFT operator '
        f'within {maximum:.3g} per transferred Fourier coordinate. Every saved event count, '
        'critical value and first-view score is independently replayed in the accompanying '
        'verification record. A threshold perturbation of plus/minus 1e-8 changes '
        f"{verification['counts_changed_by_1e8_threshold_perturbation']} count across 64 saved score arrays; "
        f"the smallest score-to-threshold distance is {verification['minimum_threshold_distance']:.3g}. "
        'Numerical agreement is not outward-rounded arithmetic.','',
        '## Interpretation','',
        'The viewing-density bound is an added modeling assumption, not a solution '
        'to the failed arbitrary-view certificate. Increasing kappa changes the '
        'critical value while the held-out alternative remains Haar; no non-Haar '
        'test sample is generated. The resulting power loss measures conservatism '
        'of this probability bound, not intrinsic sensitivity to preferred views. '
        'Power-only comparisons and all '
        'amplitude-range failures remain visible. A score fitted with knowledge of '
        'the oracle alternative does not establish practical sensitivity to unknown '
        'structural errors. Only one CTF/noise profile and the full region-removal '
        'alternative are tested; smaller deviations, unknown noise, image windows '
        'and experimentally calibrated viewing/amplitude laws remain unresolved. '
        'The construction combines classical binomial inference and nuisance-event '
        'domination. Its empirical usefulness and novelty must be assessed separately '
        'from its conditional mathematical validity. All three full Fable reviews '
        'remain rejections; no new full acceptance review occurred.',''])
    (ROOT/'research/uncertainty/paired-power-v1/MONTE-CARLO-VIEW-LAW-RESULTS.md').write_text('\n'.join(lines))
    table=[r'\begin{table*}[t]',r'\centering\footnotesize',
        r'\caption{Post-review Monte Carlo development. Projected rejection probabilities on the Haar true-map alternative against a removed-map candidate. Columns at $\kappa=1$ vary particle count; the last three use $n=10^5$. Brackets are pointwise 95\% intervals for the $\kappa=1,n=10^5$ projection. All values are rounded to three decimals; 0.000/1.000 need not be exact probabilities. The full CSV retains all 192 projections and their intervals, including correct-null controls. Increasing $\kappa$ changes the test bound, not the held-out viewing law.}',
        r'\label{tab:momentmc}',r'\begin{tabular}{lllrrrrrr}',r'\toprule',
        r'Stack & Statistic & Amplitude & $n=10^3$ & $n=10^4$ & $n=10^5$ [95\%] & $\kappa=1.01$ & $\kappa=1.1$ & $\kappa=2$ \\',r'\midrule']
    for c in r['cases']:
        if 'skipped' in c:continue
        def find(kappa,n):return next(p['true_fixed'] for p in c['projections'] if p['candidate']=='removed' and p['kappa']==kappa and p['particles']==n)
        base=[find(1,n) for n in [1000,10000,100000]];extra=[find(k,100000) for k in [1.01,1.1,2]]
        label='Power+bispectrum' if c['key'].startswith('power_bispectrum') else 'Power'
        amplitude='1' if c['amplitude_interval'][0]==c['amplitude_interval'][1] else '[.9,1.1]'
        lo,hi=base[2]['rejection_interval'];values=[f"{x['rejection_probability']:.3f}" for x in base+extra]
        values[2]+=f' [{lo:.3f}, {hi:.3f}]'
        table.append(f"{c['dataset']} & {label} & {amplitude} & "+' & '.join(values)+r' \\')
    table.extend([r'10076 & All four contrasts & Both & \multicolumn{6}{l}{Zero directions; constant scores, no simulation and no rejection.} \\',
        r'\bottomrule',r'\end{tabular}',r'\end{table*}',''])
    (ROOT/'paper/tables/moment-mc-development.tex').write_text('\n'.join(table))
    print('Wrote',len(rows),'projections and all-case report')


if __name__=='__main__':main()
