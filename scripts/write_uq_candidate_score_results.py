#!/usr/bin/env python3
"""Complete candidate-score reporting and prespecified comparison figures."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development/candidate-moment-score-v1'


def fmt(x):return '<1e-300' if x==0 else f'{x:.3g}'


def writecsv(path,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)


def main():
    records=[json.loads((BASE/ds/'summary.json').read_text()) for ds in ['10028','10049','10076']]
    baseline=json.loads((BASE/'classical-comparators.json').read_text());verified=json.loads((ROOT/'provenance/uncertainty/candidate-score-independent-verification.json').read_text())
    assert all(r['complete'] for r in records) and baseline['complete'] and verified['complete']
    rows=[];groups=[];slack=[]
    for r in records:
        for c in r['cases']:
            extra=next(e for e in baseline['cases'] if e['dataset']==r['dataset'] and e['key']==c['key'])
            slack.append(dict(dataset=r['dataset'],score=c['key'],**extra['decomposition']))
            for p in c['projections']+extra['projections']:
                for o in p['outcomes']:
                    row=dict(dataset=r['dataset'],score=c['key'],method=p['method'],kappa=p['kappa'],particles=p['particles'],deletion=o['deletion'],amplitude=o['amplitude'],
                        probability_bound=p['probability_bound'],critical_count=p['reject_at_count'],single_image_frequency=o['single_particle_probability'],rejection=o['rejection_probability'],lower=o['rejection_interval'][0],upper=o['rejection_interval'][1])
                    rows.append(row)
            for p in c['repeated_groups']+extra['repeated_groups']:
                for o in p['outcomes']:
                    groups.append(dict(dataset=r['dataset'],score=c['key'],method=p['method'],kappa=p['kappa'],particles=p['particles'],deletion=o['deletion'],amplitude=o['amplitude'],critical_count=p['reject_at_count'],groups=p['groups'],rejected=o['rejected'],fraction=o['fraction'],lower=o['pointwise_95_interval'][0],upper=o['pointwise_95_interval'][1]))
    assert len(rows)==18900 and len(groups)==6300
    writecsv(BASE/'projections.csv',rows);writecsv(BASE/'repeated-groups.csv',groups);writecsv(BASE/'variance-decomposition.csv',slack)
    def pick(ds,score,method,kappa,n,deletion,amplitude=1.):
        return next(x for x in rows if all(x[k]==v for k,v in dict(dataset=ds,score=score,method=method,kappa=kappa,particles=n,deletion=deletion,amplitude=amplitude).items()))
    lines=['# Candidate-only local-perturbation study: complete results','',
        'All three previously fitted Gaussian maps are used, with no external '
        'reference map for region or score design. Each candidate supplies a '
        '20-Angstrom smoothed peak and its own prescribed 25% local-deletion '
        'direction. Four scores compare power/combined moments and raw/scale-'
        'orthogonal directions. Training uses 4,096 fresh Haar views, calibration '
        '32,768 views with two groups of 32 noises, and testing 131,072 independent '
        'Haar views per stack. All noise, transfer and viewing assumptions remain '
        'specified simulators. These are new simulations of fitted candidates, '
        'not new reconstructions or applications to experimental test particles.','',
        '## Paired-variance projections at n=10,000 and amplitude one','',
        'Each entry is a binomial projection, not observed repeated-dataset power. '
        'The complete CSV includes pointwise 95% intervals, both other amplitudes, '
        'every sample size and all seven calibration procedures. Rounded tiny values '
        'are displayed in scientific notation. Exact non-rejection and numerical '
        'underflow are distinguished by the saved critical counts.','',
        '| Stack | Score | kappa | 10% removal | 25% | 50% | 100% |',
        '|---|---|---:|---:|---:|---:|---:|']
    for r in records:
        ds=r['dataset']
        for c in r['cases']:
            for kappa in [1.,1.1]:
                values=[pick(ds,c['key'],'view_variance',kappa,10000,d)['rejection'] for d in [.1,.25,.5,1.]]
                lines.append(f'| {ds} | {c["key"]} | {kappa:g} | '+' | '.join(map(fmt,values))+' |')
    lines+=['','## Classical comparison on full deletion, kappa=1.1, n=10,000','',
        'These separately declared tests share draws. No minimum across methods '
        'is used as a valid combined test. The table is restricted by score type '
        '(the scale-orthogonal arm); all raw-arm results remain in the CSV.','',
        '| Stack | Score | Paired variance | CVaR split | CVaR DKW | Grouped ratio | Unpaired variance |',
        '|---|---|---:|---:|---:|---:|---:|']
    for ds in ['10028','10049','10076']:
        for score in ['power_scale_orthogonal','power_bispectrum_scale_orthogonal']:
            values=[pick(ds,score,m,1.1,10000,1.)['rejection'] for m in ['view_variance','cvar_split','cvar_dkw','grouped_ratio','unpaired_variance']]
            lines.append(f'| {ds} | {score} | '+' | '.join(map(fmt,values))+' |')
    nulls=[r for r in rows if r['deletion']==0 and r['particles']==10000]
    null_groups=[r for r in groups if r['deletion']==0]
    maxnull=max(nulls,key=lambda x:x['rejection']);maxobserved=max(null_groups,key=lambda x:x['rejected'])
    lines+=['','## Correct-null controls, scope and verification','',
        f'Across every method/score/kappa/amplitude at n=10,000, the largest correct-'
        f'null point projection is {maxnull["rejection"]:.6g}. Across all '
        f'{len(null_groups)} correct-null repeated-group cells, the largest observed '
        f'count is {maxobserved["rejected"]}/128 (pointwise exact 95% interval '
        f'[{maxobserved["lower"]:.5g}, {maxobserved["upper"]:.5g}]). '
        'Calibration is fixed, cells share draws, and these cells must not be '
        'pooled as independent repetitions. This is not a calibration-repeat '
        'or unconditional-error experiment.','',
        f'The full report retains {len(rows):,} scalar projections (1,260 method rows) '
        f'and {len(groups):,} actual group-outcome cells. Both smaller deletions and '
        'all zero-deletion/amplitude controls remain visible. The removed-region '
        'L2 fractions are '+', '.join(f'{r["dataset"]}: {r["region"]["removed_L2_fraction"]:.5f}' for r in records)+
        '; multiplying these by the deletion amount gives the relative L2 change '
        'in the candidate. Normalization and the transfer/noise fixture are '
        'specified, not estimated experimental signal-to-noise ratios.','',
        'Independent FFT convolution reproduces all three candidate-only region '
        'centers. Separate feature arithmetic and least-squares projection reproduce '
        f'all 12 directions within {verified["maximum_training_direction_difference"]:.3g}. '
        f'The check includes {verified["calibration_first_polynomials"]} first-view '
        'calibration cubics, 180 direct held-out scores, every held-out event count, '
        f'{verified["critical_value_checks"]} minimal critical values, all 6,300 group '
        'vectors and all 18,900 rejection projections. Maximum direct-score '
        f'discrepancy is {verified["maximum_direct_score_difference"]:.3g}. '
        'It is a numerical replay, not a global floating-point certificate.','',
        'The three Mac processes take '+', '.join(f'{r["dataset"]}: {r["seconds"]:.2f}s' for r in records)+
        '. They overlap; summed process durations are not elapsed wall time.','',
        '## Scientific interpretation','',
        'Candidate-driven design removes the need to know an external true-map '
        'alternative for training this score. It does not remove the hand-specified '
        'deletion family or calibrate experimental imaging nuisance. Scale '
        'orthogonalization improves several full/half-deletion cases, while '
        '10--25% changes remain weak. Stronger view-law protection again removes '
        'much of the sensitivity. Unlike the previous selected oracle contrasts, '
        'these candidate scores show substantial empirical between-view variation; '
        'the paired variance decomposition is retained for every score. On 10028 '
        'and 10049 the CVaR comparator can be tighter than the paired variance '
        'method. This reversal is retained and rules out a general empirical '
        'superiority claim from the earlier experiment.','',
        'A rejection concerns the complete candidate under the stated simulator. '
        'Other structural/imaging errors can change the same moment score, so '
        'neither rejection nor non-rejection is a unique local-error attribution '
        'or confidence interval for regional occupancy. All three full Fable '
        'reviews still reject. No new full acceptance assessment has occurred.','']
    (ROOT/'research/uncertainty/paired-power-v1/CANDIDATE-SCORE-RESULTS.md').write_text('\n'.join(lines))
    # All three stacks; scale-orthogonal power and combined moments; no outcome selection.
    fig,axes=plt.subplots(2,3,figsize=(10.2,5.9),sharex=True,sharey=True,constrained_layout=True)
    for col,ds in enumerate(['10028','10049','10076']):
        for row,n in enumerate([10000,100000]):
            ax=axes[row,col]
            for family,color in [('power','#2166ac'),('power_bispectrum','#b35806')]:
                score=family+'_scale_orthogonal'
                for kappa,style in [(1.,'-'),(1.1,'--')]:
                    selected=[pick(ds,score,'view_variance',kappa,n,d) for d in [0,.1,.25,.5,1.]]
                    x=np.array([p['deletion'] for p in selected]);y=np.array([p['rejection'] for p in selected]);lo=np.array([p['lower'] for p in selected]);hi=np.array([p['upper'] for p in selected])
                    label=('Power' if family=='power' else 'Power + bispectrum')+f', k={kappa:g}'
                    ax.plot(x,y,style,color=color,marker='o',markersize=3,lw=1.3,label=label);ax.fill_between(x,lo,hi,color=color,alpha=.07)
            ax.set_title(f'EMPIAR {ds}, n={n:,}',fontsize=10);ax.set_ylim(-.02,1.03);ax.set_xticks([0,.25,.5,1.]);ax.grid(alpha=.18)
            if col==0:ax.set_ylabel('Projected rejection probability')
            if row==1:ax.set_xlabel('Fraction of nominated region removed')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=4,fontsize=8)
    fig.suptitle('Candidate-derived scores: Haar test views, amplitude 1\nPointwise intervals; calibration uncertainty not included',fontsize=11)
    dest=ROOT/'paper/figures/candidate-score-projections.pdf';fig.savefig(dest);fig.savefig(dest.with_suffix('.png'),dpi=160);plt.close(fig)
    print('Wrote',len(rows),'projections,',len(groups),'group outcomes and figure')


if __name__=='__main__':main()
