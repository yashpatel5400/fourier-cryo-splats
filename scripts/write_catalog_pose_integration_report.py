#!/usr/bin/env python3
"""All cases from the final allowed integration repair, including failures."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def main():
    out=BASE/'catalog-pose-integration-report-v1'
    if out.exists():raise ValueError('Preserve every report attempt')
    entries=[]
    autopsies={name:json.loads((BASE/f'{name}-pose-autopsy-v1/summary.json').read_text()) for name in ['adaptive','catalog']}
    assert all(v['complete'] for v in autopsies.values())
    for ds in ['10028','10049','10076']:
        d=BASE/'catalog-pose-integration-v1'/ds
        j=json.loads((d/'summary.json').read_text());v=json.loads((d/'independent-check.json').read_text())
        old=json.loads((BASE/'adaptive-pose-integration-v1'/ds/'summary.json').read_text())
        assert j['complete'] and v['complete'] and old['complete'];entries.append((ds,j,v,old))
    out.mkdir();rows=[];table=[];errors=[];checks=[];autopsy_table=[];brackets=[]
    fig,axes=plt.subplots(2,3,figsize=(11,6.5),constrained_layout=True)
    for col,(ds,j,v,old) in enumerate(entries):
        dif=[];ess=[];integral_difference=[];prefix=[]
        for c,o in zip(j['cases'],old['cases']):
            assert (c['position'],c['state'],c['index'])==(o['position'],o['state'],o['index'])
            bank=[b['levels'][-1] for b in c['banks']]
            dif.append(c['between_bank_logratio_difference']);ess.append([b['ess'] for b in bank])
            integral_difference.append(np.asarray(bank[0]['log_integrals'])-bank[1]['log_integrals'])
            prefix.append([b['levels'][-1]['log_ratio']-b['levels'][0]['log_ratio'] for b in c['banks']])
            row=dict(dataset=ds,position=c['position'],state=c['state'],source_index=c['index'],
                     old_bank_difference=o['between_bank_logratio_difference'],new_bank_difference=dif[-1])
            for b in [0,1]:
                row[f'bank{b}_logratio']=bank[b]['log_ratio']
                for model in [0,1]:
                    row[f'bank{b}_model{model}_ess']=bank[b]['ess'][model]
                    row[f'bank{b}_model{model}_log_integral']=bank[b]['log_integrals'][model]
            rows.append(row)
        dif=np.array(dif);ess=np.array(ess);absolute=np.array(integral_difference);prefix=np.array(prefix)
        med=np.median(ess,axis=0);low=np.quantile(ess,.1,axis=0)
        table.append(f"| {ds} | {np.sum(dif<=.01)}/128 ({np.mean(dif<=.01):.4f}) | {med.min():.1f}–{med.max():.1f} | {low.min():.1f}–{low.max():.1f} | {'pass' if j['gate']['passed'] else 'fail'} | {j['seconds']/60:.2f} |")
        errors.append(f"| {ds} | {np.sqrt(np.mean(dif**2)):.6g} | {dif.max():.6g} | {np.sqrt(np.mean(absolute**2)):.6g} | {np.sqrt(np.mean(prefix**2)):.6g} |")
        olddiff=np.array([r['between_bank_logratio_difference'] for r in old['cases']])
        ax=axes[0,col];ax.plot(np.arange(128),olddiff,'.',color='#999999',alpha=.6,ms=3,label='First proposal')
        ax.plot(np.arange(128),dif,'.',color='#2166ac',ms=4,label='Final repair');ax.axhline(.01,color='#b2182b',ls='--',lw=1)
        ax.set(yscale='log',xlabel='Image case (all 128)',ylabel='Absolute bank log-ratio difference',title=f'EMPIAR {ds}');ax.grid(alpha=.15)
        for bank,style in [(0,'-'),(1,'--')]:
            for model,color in [(0,'#2166ac'),(1,'#b2182b')]:
                axes[1,col].plot(np.arange(1,129)/128,np.sort(ess[:,bank,model]),ls=style,color=color,label=f'Bank {bank}, model {model}')
        axes[1,col].axhline(256,color='black',lw=1);axes[1,col].set(yscale='log',xlabel='Empirical quantile',ylabel='Importance ESS');axes[1,col].grid(alpha=.15)
        if col==0:axes[0,col].legend(fontsize=8)
        checks.append(f"- {ds}: all 128 integral/ESS/gate summaries and observations replay. Independent full-4D covariance densities check {v['checked_densities']:,} points (32 fixed points in each bank per image), maximum log-density error {max(r['maximum_log_density_error'] for r in v['images']):.3g}. Six direct 64³-cell residual calculations cover both models at all 220 frequencies, maximum discrepancy {max(r['error'] for r in v['physical_checks']):.3g}.")
        for attempt in ['adaptive','catalog']:
            auto=next(d for d in autopsies[attempt]['datasets'] if d['dataset']==ds)
            for s in auto['failure_strata']:
                stats=s['bank_statistics']
                autopsy_table.append(f"| {ds} | {attempt} | {'within' if s['ratio_within_tolerance'] else 'outside'} | {s['images']} | {stats['haar_weight_mass']['mean']:.6g} | {stats['planted_jump']['quantiles']['1']:.6g} | {s['median_standardized_discrepancy']:.4g} |")
            assert all(x['clipped_eigenvalues']==0 for x in auto['images'])
            widths=[x['width'] for x in auto['stack_brackets']]
            brackets.append(f"| {ds} | {attempt} | {min(widths):.6g}–{max(widths):.6g} |")
    fig.legend(*axes[1,0].get_legend_handles_labels(),loc='outside lower center',ncol=4,fontsize=9)
    for ext in ['png','pdf']:fig.savefig(out/f'final-integration-gate.{ext}',dpi=180)
    plt.close(fig)
    with (out/'all-images.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    overall=all(j['gate']['passed'] for _,j,_,_ in entries)
    manifest=dict(complete=True,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases=len(rows),gate_passed=overall)
    (out/'summary.json').write_text(json.dumps(manifest,indent=2)+'\n')
    lines=['# Final catalogue repair of orientation integration','',
        '1 October 2026 UTC. **All three stacks pass the unchanged numerical gate.** This completes the single permitted substantive repair. All 384 image cases, both banks, original optimizer failures and thirteen remaining out-of-tolerance images are retained. The cohort, maps, 8,192 draws per bank, tolerance and median-ESS thresholds were fixed before outcomes; the catalogue proposal uses continuous rotations and a full normalized mixture density.',
        '', '| EMPIAR | Images within .01 | Four median ESS range | Four lower-decile ESS range | Stack gate | Minutes |',
        '|---|---:|---:|---:|---|---:|',*table,'',
        'Each stack must have at least 90% of its images within .01 log-ratio units between banks and each of its four median ESS values at least 256. Lower-decile ESS is secondary: the later consultation suggested 900/450 on 10049/10076, but these were not substituted into the frozen gate. Timings are concurrent Mac wall times, not a GPU or matched-baseline speedup. No external compute was rented.',
        '', '| EMPIAR | RMS ratio difference | Maximum ratio difference | RMS absolute log-integral difference | RMS prefix ratio change |',
        '|---|---:|---:|---:|---:|',*errors,'',
        'Two proposals can miss the same mass; ratio error also benefits from common-mode cancellation between the two nearby maps. Passing this per-image gate does not certify all parameter values, a product likelihood over a large dataset, or real-image marginal likelihoods. It does not validate Haar experimental views, white noise, known amplitudes, correct maps or uncertainty coverage. The [population materiality screen](POPULATION-MATERIALITY-RESULTS.md) fails independently, so this numerical pass does not authorize the proposed scientific method.',
        '', '## Independent checks','',*checks,'',
        'The replay recomputes every stored importance integral and verifies identical observations to attempt one. The density check uses explicitly transformed inverse 4D covariance matrices and the physical check explicitly sums cells; neither calls the runner’s corresponding shortcut. Only the stated density/physical subset is checked independently, not every integrand point.',
        '', '## Post hoc saved-draw autopsy','',
        '| EMPIAR | Proposal | Original .01 criterion | Images | Mean matching-state Haar weight mass | Largest planted log-integral jump | Median bank discrepancy / delta SE |',
        '|---|---|---|---:|---:|---:|---:|',*autopsy_table,'',
        'The planted diagnostic replaces draw zero with the saved generating pose under the matching state only. Its [reciprocal-unbiased identity and assumptions](PLANTED-TRUTH-INTEGRATION-AUDIT.md) are proved by exchangeability and checked by exact finite-space enumeration. Large jumps expose sensitivity to posterior mass absent from ordinary saved draws; they do not reveal the exact integral. In particular, the old 10076 within-tolerance group contains a jump of 1.50 log units, illustrating why bank agreement is insufficient. After repair its largest within-tolerance jump is .0626. No original mode has a clipped Hessian eigenvalue, so clipping does not explain these failures. Delta SEs remain tail-sensitive plug-in diagnostics.',
        '', '| EMPIAR | Proposal | Four 64-image stochastic bracket widths |',
        '|---|---|---:|',*brackets,'',
        'These are pointwise 95% stochastic brackets from one forward and one planted reciprocal estimate per independent image, with .025 allocated to each tail. Each bracket concerns one generating state and one bank. The paired states are never multiplied together, and no simultaneous coverage across the displayed brackets is asserted. The additive floor 2 log(40) is about 7.378 log units. These intervals are simulation-only and too coarse to act as a general numerical certificate.',
        '', '## Decision and chronology','',
        'The numerical branch ends here as planned: no further proposal, images, draws or favorable subset are selected. The method consultation was complete at the provider at 11:38:56 UTC but had not been retrieved or read when the 11:40 scheduling amendment and run were made. Its later recommendation of a local-majority proposal differs from the frozen .45 local/.50 catalogue/.05 Haar choice. The [response correction](reviews/post-round04-method-consultation/response.md) preserves that distinction. Original protocols and provider text remain unchanged.',
        '', 'All numerical points, labels, catalogue weights, physical kernels, source hashes and unfavorable cases are preserved. This is a classical integration repair and a development audit, not a new uncertainty theorem, an experimental population estimate or a fifth full-paper review.']
    (ROOT/'research/uncertainty/CATALOG-POSE-INTEGRATION-RESULTS.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main()
