#!/usr/bin/env python3
"""Complete tables and figures from the declared post hoc saved-array replay."""
import csv
import hashlib
import json
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fourier_splats.uq_end_to_end import binomial_interval

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development/refitting-bias-reanalysis-v1'
OUT=ROOT/'results/uncertainty/development/refitting-bias-summary-v1'
KEYS=['dataset','template','target','method','image_mode']


def read_rows(path):
    with path.open() as f: rows=list(csv.DictReader(f))
    for r in rows:
        for k,v in r.items():
            if k in ['dataset','template','target','method','image_mode']:continue
            r[k]= v=='True' if v in ['True','False'] else float(v)
    return rows


def stats(a):
    a=np.asarray(a,float)
    return dict(mean=float(a.mean()),se=float(a.std(ddof=1)/np.sqrt(len(a))),
        median=float(np.median(a)),minimum=float(a.min()),maximum=float(a.max()))


def main():
    OUT.mkdir(parents=True,exist_ok=True); groups=defaultdict(list); diags=defaultdict(list); metas={}
    for ds in ['10028','10049','10076']:
        p=BASE/ds; meta=json.loads((p/'metadata.json').read_text());metas[ds]=meta
        assert meta['complete'] and meta['replicates']==200 and meta['rows']==10400 and meta['realized_envelopes']==800
        for name,h in meta['outputs'].items():assert hashlib.sha256((p/name).read_bytes()).hexdigest()==h
        for r in read_rows(p/'rows.csv'):groups[tuple(r[k] for k in KEYS)].append(r)
        for r in read_rows(p/'realized-envelopes.csv'):diags[tuple(r[k] for k in ['dataset','template','target'])].append(r)
    summaries=[]; lookup={}
    for key,rows in groups.items():
        assert len(rows)==200 and len(set(r['replicate'] for r in rows))==200
        d=dict(zip(KEYS,key));d['n']=200;d['truth']=rows[0]['truth'];d['pilot']=rows[0]['pilot'];d['pilot_error']=rows[0]['pilot_error']
        d['rmse']=float(np.sqrt(np.mean([r['error']**2 for r in rows])))
        d['pilot_absolute_error']=abs(d['pilot_error']);d['rmse_over_pilot_error']=d['rmse']/d['pilot_absolute_error']
        for metric in ['error','standardized_error','realized_signal_bias','standardized_signal_bias','standardized_projected_noise','noise_sd','original_half_width']:
            for name,value in stats([r[metric] for r in rows]).items():d[metric+'_'+name]=value
        for tag in ['original','B0','B0.5','B1','B2']:
            successes=sum(r[tag+'_covered'] for r in rows);ci=binomial_interval(successes,200)
            d[tag+'_coverage']=successes/200;d[tag+'_coverage_count']=successes;d[tag+'_coverage_lower']=ci[0];d[tag+'_coverage_upper']=ci[1]
            d[tag+'_half_width_median']=float(np.median([r[tag+'_half_width'] for r in rows]))
            if tag!='original':d[tag+'_contains_generator']=rows[0][tag+'_contains_generator']
        summaries.append(d);lookup[key]=d
    pairs=[]
    for key,rows in groups.items():
        if key[-1]!='same_image':continue
        other={r['replicate']:r for r in groups[key[:-1]+('independent_image',)]}
        differences=np.array([r['center']-other[r['replicate']]['center'] for r in rows])
        z=np.array([(r['center']-other[r['replicate']]['center'])/r['noise_sd'] for r in rows])
        p=dict(zip(KEYS[:-1],key[:-1]));p['n']=200
        for label,values in [('paired_center_difference',differences),('paired_standardized_difference',z)]:
            p.update({label+'_'+k:v for k,v in stats(values).items()})
        pairs.append(p)
    envelopes=[]
    for key,rows in diags.items():
        d=dict(zip(['dataset','template','target'],key));d['n']=len(rows)
        d['inside_ball_count']=sum(r['pose_inside_calibrated_ball'] for r in rows)
        for tag in ['realized_total_bias_envelope','realized_pose_bias_envelope','original_deterministic_bias','original_mixed_half_width','first_order_only_half_width','first_order_tail','fixed_bias']:
            values=[r[tag]/r['no_data_half_width'] for r in rows]
            d.update({tag+'_relative_'+k:v for k,v in stats(values).items()})
        d['first_order_without_remainder_fallback_count']=sum(r['first_order_only_half_width']>r['no_data_half_width'] for r in rows)
        envelopes.append(d)
    for name,rows in [('groups',summaries),('paired-differences',pairs),('realized-envelopes',envelopes)]:
        with (OUT/(name+'.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=dict(complete=True,scope='Post hoc reanalysis of saved data; no new simulated observations and no prospective comparison.',
        metadata=metas,groups=summaries,paired_differences=pairs,realized_envelopes=envelopes,
        counts=dict(datasets=600,estimator_image_rows=sum(len(g) for g in groups.values()),groups=len(groups),realized_envelopes=sum(len(g) for g in diags.values())))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,3,figsize=(7.1,4.05),sharey='row')
    templates=['true_pose','oracle_reference','independent_pilot'];positions=np.arange(3)
    for col,ds in enumerate(['10028','10049','10076']):
        for target,color in [('center','#0072B2'),('contrast','#D55E00')]:
            for mode,marker,offset in [('same_image','o',-.08),('independent_image','s',.08)]:
                rows=[lookup[ds,t,target,'fixed_folded',mode] for t in templates]
                extra=-.09 if target=='center' else .09
                label=('Average' if target=='center' else 'Contrast')+(' · same image' if mode=='same_image' else ' · independent image')
                axes[0,col].errorbar(positions+extra+offset*.45,[r['standardized_error_mean'] for r in rows],
                    yerr=[1.96*r['standardized_error_se'] for r in rows],color=color,marker=marker,
                    mfc=color if mode=='same_image' else 'white',ms=3.5,lw=.8,linestyle='-' if mode=='same_image' else '--',label=label)
                cov=np.array([r['B0_coverage'] for r in rows]);lower=np.array([r['B0_coverage_lower'] for r in rows]);upper=np.array([r['B0_coverage_upper'] for r in rows])
                axes[1,col].errorbar(positions+extra+offset*.45,cov,yerr=[cov-lower,upper-cov],color=color,marker=marker,
                    mfc=color if mode=='same_image' else 'white',ms=3.5,lw=.8,linestyle='-' if mode=='same_image' else '--')
        axes[0,col].axhline(0,color='.6',lw=.7);axes[0,col].set_title('EMPIAR-'+ds)
        axes[1,col].axhline(.95,color='.4',ls=':',lw=.8)
        for ax in axes[:,col]:ax.set_xticks(positions,['True\nposes','Oracle\ntemplate','Pilot\ntemplate']);ax.grid(axis='y',alpha=.15)
    axes[0,0].set_ylabel('Mean error / noise SD');axes[1,0].set_ylabel('Noise-only coverage')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,-.008))
    fig.tight_layout(rect=[0,.095,1,1]);figdir=ROOT/'paper/figures'
    fig.savefig(figdir/'refitting-alignment-bias.pdf',bbox_inches='tight');fig.savefig(figdir/'refitting-alignment-bias.png',dpi=180,bbox_inches='tight');plt.close(fig)
    report=['# Saved-array alignment diagnosis after review 3','',
        'Post hoc analysis of all 600 frozen trials. All numerical inputs are hash-checked; no new random observations are generated. The true-pose controls are newly fitted deterministic operators on the saved images. The reduced-radius rows retain their original weights and explicitly flag whether their class contains the generator. All individual Monte Carlo intervals are nonsimultaneous.','',
        '## Audit-weight results','',
        '| Stack | Alignment | Target | Mean error / SD, same / independent | Noise-only coverage, same / independent | RMSE / pilot error, same / independent |',
        '|---|---|---|---|---|---|']
    for ds in ['10028','10049','10076']:
        for t in templates:
            for target in ['center','contrast']:
                a=lookup[ds,t,target,'fixed_folded','same_image'];b=lookup[ds,t,target,'fixed_folded','independent_image']
                report.append(f"| {ds} | {t} | {target} | {a['standardized_error_mean']:.3f} / {b['standardized_error_mean']:.3f} | {a['B0_coverage']:.3f} / {b['B0_coverage']:.3f} | {a['rmse_over_pilot_error']:.3f} / {b['rmse_over_pilot_error']:.3f} |")
    report.extend(['','## Continuous class at realized poses','','These are oracle diagnostics at the realized pose pair, not deployable intervals or uniform pose-set bounds. All norms use the continuous cube; Gauss quadrature errors are bounded analytically, while NUFFT and rounding are only checked numerically.','','| Stack | Template | Target | Total bias / no data, median | Pose-only bias / no data, median | First-order half-width / no data, median |','|---|---|---|---|---|---|'])
    for d in envelopes:
        report.append(f"| {d['dataset']} | {d['template']} | {d['target']} | {d['realized_total_bias_envelope_relative_median']:.3f} | {d['realized_pose_bias_envelope_relative_median']:.3f} | {d['first_order_only_half_width_relative_median']:.3f} |")
    report.extend(['','## Exact pilot-to-generator distance',''])
    for ds,m in metas.items():report.append(f"- {ds}: {m['density_distance']:.6f}; all 6 true-pose fits converged: {m['true_pose_fits_converged']}; maximum replay discrepancy {m['max_saved_center_discrepancy']:.3g}.")
    report.extend(['','The machine-readable summaries retain all 156 estimator/image cells, 78 paired comparisons and 2,400 realized envelope calculations. Class-radius and noise-only coverage are diagnostic, not evidence that these smaller classes are valid experimentally.',''])
    (ROOT/'research/uncertainty/REFITTING-BIAS-REANALYSIS-RESULTS.md').write_text('\n'.join(report))
    print(json.dumps(summary['counts']))

if __name__=='__main__':main()
