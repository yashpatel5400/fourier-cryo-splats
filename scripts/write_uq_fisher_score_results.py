#!/usr/bin/env python3
"""All outcomes and conservative pointwise comparisons for the fresh Fisher study."""
import csv,json
from pathlib import Path
import numpy as np
from scipy.stats import beta,binom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development/candidate-fisher-score-v1'
def writecsv(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def fmt(x):return '<1e-300' if x==0 else f'{x:.3g}'
def main():
    records=[json.loads((BASE/ds/'summary.json').read_text()) for ds in ['10028','10049','10076']]
    assert all(r['complete'] for r in records)
    verified=json.loads((ROOT/'provenance/uncertainty/fisher-score-independent-verification.json').read_text());assert verified['complete']
    rows=[];groups=[];pairs=[];designs=[]
    for r in records:
        ds=r['dataset']
        with np.load(ROOT/r['array']['path']) as f:
            scores=f['heldout_scores'];cov=f['training_covariance'];nq=len(f['q']);events=[]
            for j,c in enumerate(r['cases']):
                w=f[c['key']+'_direction'];p=len(w);reg=.9*cov[:p,:p]+.1*np.trace(cov[:p,:p])/p*np.eye(p)
                gap=c['design']['training_design_alternative_mean']-c['design']['training_null_mean'];sv=scores[j,0,1].var(ddof=1)
                testgap=scores[j,2,1].mean()-scores[j,0,1].mean()
                designs.append(dict(dataset=ds,score=c['key'],training_gap=gap,regularized_training_variance=w@reg@w,
                    training_squared_separation=gap**2/(w@reg@w),test_null_mean=scores[j,0,1].mean(),
                    test_design_gap=testgap,test_null_variance=sv,test_squared_separation=testgap**2/sv))
                events.append(scores[j]>c['design']['threshold'])
                for p in c['projections']:
                    for o in p['outcomes']:
                        rows.append(dict(dataset=ds,score=c['key'],method=p['method'],kappa=p['kappa'],particles=p['particles'],
                            deletion=o['deletion'],amplitude=o['amplitude'],probability_bound=p['probability_bound'],
                            critical_count=p['reject_at_count'],single_image_frequency=o['single_particle_probability'],
                            rejection=o['rejection_probability'],lower=o['rejection_interval'][0],upper=o['rejection_interval'][1]))
                for p in c['repeated_groups']:
                    for o in p['outcomes']:
                        groups.append(dict(dataset=ds,score=c['key'],method=p['method'],kappa=p['kappa'],particles=p['particles'],
                            deletion=o['deletion'],amplitude=o['amplitude'],critical_count=p['reject_at_count'],groups=p['groups'],
                            rejected=o['rejected'],fraction=o['fraction'],lower=o['pointwise_95_interval'][0],upper=o['pointwise_95_interval'][1]))
            for family,j in [('power',0),('power_bispectrum',2)]:
                mc,fc=r['cases'][j:j+2];assert mc['key']==family+'_matched' and fc['key']==family+'_fisher'
                for pm,pf in zip(mc['projections'],fc['projections']):
                    assert all(pm[k]==pf[k] for k in ['method','kappa','particles'])
                    for om,of in zip(pm['outcomes'],pf['outcomes']):
                        di=r['deletions'].index(om['deletion']);ai=r['amplitudes'].index(om['amplitude']);em,ef=events[j][di,ai],events[j+1][di,ai]
                        cm,cf=int(em.sum()),int(ef.sum());total=r['heldout_views'];n=pm['particles'];km,kf=pm['reject_at_count'],pf['reject_at_count']
                        def ci(c):return (0. if c==0 else beta.ppf(.0125,c,total-c+1),1. if c==total else beta.ppf(.9875,c+1,total-c))
                        lm,um=ci(cm);lf,uf=ci(cf)
                        pairs.append(dict(dataset=ds,family=family,method=pm['method'],kappa=pm['kappa'],particles=n,
                            deletion=om['deletion'],amplitude=om['amplitude'],matched=om['rejection_probability'],fisher=of['rejection_probability'],
                            difference=of['rejection_probability']-om['rejection_probability'],
                            lower=binom.sf(kf-1,n,lf)-binom.sf(km-1,n,um),upper=binom.sf(kf-1,n,uf)-binom.sf(km-1,n,lm),
                            neither=int(np.sum(~em&~ef)),matched_only=int(np.sum(em&~ef)),fisher_only=int(np.sum(~em&ef)),both=int(np.sum(em&ef))))
        print(ds,'all paired comparison records complete',flush=True)
    assert len(rows)==18900 and len(groups)==6300 and len(pairs)==9450
    for name,values in [('projections',rows),('repeated-groups',groups),('fisher-minus-matched',pairs),('score-design',designs)]:writecsv(BASE/(name+'.csv'),values)
    def pick(ds,fam,method,kappa,n,deletion):return next(r for r in pairs if all(r[k]==v for k,v in dict(dataset=ds,family=fam,method=method,kappa=kappa,particles=n,deletion=deletion,amplitude=1.).items()))
    lines=['# Fresh constrained Fisher comparison: complete results','',
        'All three candidates use 65,536 new training views, 32,768 new calibration views '
        'with two groups of 32 noises, and 131,072 independent Haar test views. '
        'Matched and constrained Fisher scores share the same training mean gap; '
        'Fisher uses covariance shrinkage 0.1 without tuning. All twelve scores '
        'are retained. Seven calibration methods, five viewing caps, three particle '
        'counts, five deletions and three amplitudes are in the CSVs.','',
        '## Detection at n=10,000, amplitude one','',
        'Entries are matched / Fisher rejection projections. This table includes both '
        'paired variance and split CVaR; neither score nor calibration method is '
        'selected by an unadjusted minimum. Other methods, amplitudes, caps and '
        'particle counts remain in the full archive.','',
        '| Stack | Family | Bound | kappa | 10% deletion | 25% | 50% | 100% |','|---|---|---|---:|---:|---:|---:|---:|']
    for ds in ['10028','10049','10076']:
        for fam in ['power','power_bispectrum']:
            for method in ['view_variance','cvar_split']:
                for kappa in [1.,1.1]:
                    cell=[pick(ds,fam,method,kappa,10000,d) for d in [.1,.25,.5,1.]]
                    lines.append(f'| {ds} | {fam} | {method} | {kappa:g} | '+' | '.join(fmt(p['matched'])+' / '+fmt(p['fisher']) for p in cell)+' |')
    lines+=['','## Pointwise difference intervals','',
        'Each row below compares 25% deletion at n=100,000, amplitude one, kappa=1. '
        'Intervals combine two exact 97.5% single-image intervals by a union bound '
        'and monotone binomial transforms. They have pointwise at least 95% '
        'coverage conditional on frozen training/calibration under the simulator. '
        'They do not quantify training/calibration variability and are not simultaneous. '
        'All 9,450 pairs and their joint event tables are retained.','',
        '| Stack | Family | Method | Matched | Fisher | Difference | Interval |','|---|---|---|---:|---:|---:|---|']
    for ds in ['10028','10049','10076']:
        for fam in ['power','power_bispectrum']:
            for method in ['view_variance','cvar_split']:
                p=pick(ds,fam,method,1.,100000,.25)
                lines.append(f'| {ds} | {fam} | {method} | {fmt(p["matched"])} | {fmt(p["fisher"])} | {p["difference"]:.3g} | [{p["lower"]:.3g}, {p["upper"]:.3g}] |')
    nulls=[r for r in rows if r['deletion']==0 and r['particles']==10000];worst=max(nulls,key=lambda r:r['rejection'])
    controls=[r for r in groups if r['deletion']==0];g=max(controls,key=lambda r:r['rejected'])
    lines+=['','## Correct-null controls and verification','',
        f'Largest n=10,000 correct-null projection: {worst["rejection"]:.8g}, interval [{worst["lower"]:.8g}, {worst["upper"]:.8g}], at {worst["dataset"]}, {worst["score"]}, {worst["method"]}, kappa={worst["kappa"]}, amplitude={worst["amplitude"]}.',
        f'Largest correct-null actual group count: {g["rejected"]}/128, pointwise interval [{g["lower"]:.5g}, {g["upper"]:.5g}]. '
        'Calibration is fixed, cells share draws, and the 1,260 null cells cannot be pooled as independent repetitions.','',
        f'Independent replay checks {verified["covariance_checks"]} full feature covariances, '
        f'{verified["training_direction_checks"]} directions via a null-space optimization, '
        f'{verified["bound_checks"]} probability bounds, {verified["critical_value_checks"]} minimal critical values, '
        f'{verified["projection_checks"]} projections and {verified["group_vector_checks"]} group vectors. '
        'Physical first views, first calibration cubics and all training noise draws are also checked. '
        'This does not certify all physical operator calls or rounding.','',
        'Whole-candidate testing still depends on a hand-specified deletion family and '
        'known transfer/noise/viewing assumptions. It does not localize a rejection '
        'to the nominated region or prove experimental regional-density coverage. '
        'This classical baseline does not supply a new optimization theorem or an ICML acceptance assessment.','']
    (ROOT/'research/uncertainty/paired-power-v1/FISHER-SCORE-RESULTS.md').write_text('\n'.join(lines))
    plt.rcParams.update({'font.size':9});fig,axes=plt.subplots(2,3,figsize=(10.5,6.3),sharex=True,sharey=True)
    for col,ds in enumerate(['10028','10049','10076']):
        for row,method in enumerate(['view_variance','cvar_split']):
            ax=axes[row,col]
            for family,color in [('power','#1261a0'),('power_bispectrum','#c15c0b')]:
                for variant,style in [('matched','--'),('fisher','-')]:
                    a=[pick(ds,family,method,1.1,10000,d) for d in [0.,.1,.25,.5,1.]]
                    # Show point estimates; complete individual and difference intervals in CSV.
                    ax.plot([p['deletion'] for p in a],[p[variant] for p in a],style,marker='o',ms=3,color=color,label=family.replace('_',' + ')+' '+variant)
            ax.set_title(f'{ds}, '+('paired variance' if row==0 else 'split CVaR'));ax.grid(alpha=.2);ax.set_ylim(-.02,1.02)
            if col==0:ax.set_ylabel('Projected rejection probability')
            if row==1:ax.set_xlabel('Fraction of nominated region removed')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,bbox_to_anchor=(.5,-.025))
    fig.suptitle('Fresh classical score comparison: n=10,000, amplitude 1, viewing cap 1.1\nConditional simulator projections; intervals and all cases retained in CSV')
    fig.tight_layout(rect=[0,.065,1,.94]);out=ROOT/'paper/figures/fisher-score-projections'
    fig.savefig(out.with_suffix('.pdf'),bbox_inches='tight');fig.savefig(out.with_suffix('.png'),dpi=150,bbox_inches='tight')
    print('COMPLETE',len(rows),len(groups),len(pairs),flush=True)
if __name__=='__main__':main()
