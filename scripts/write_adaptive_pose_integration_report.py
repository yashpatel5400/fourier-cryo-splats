#!/usr/bin/env python3
"""Report the prespecified gate on every image, including optimizer failures."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def main():
    out=BASE/'adaptive-pose-integration-report-v1'
    if out.exists():raise ValueError('Preserve report attempts')
    datasets=[]
    for ds in ['10028','10049','10076']:
        j=json.loads((BASE/'adaptive-pose-integration-v1'/ds/'summary.json').read_text())
        v=json.loads((BASE/'adaptive-pose-integration-v1'/ds/'independent-check.json').read_text())
        assert j['complete'] and v['complete'];datasets.append((ds,j,v))
    out.mkdir();rows=[];table=[];fitrows=[];secondary=[]
    fig,axes=plt.subplots(2,3,figsize=(11,6.2),constrained_layout=True)
    for col,(ds,j,v) in enumerate(datasets):
        dif=[];ess=[];density_diff=[];prefix=[]
        with np.load(BASE/'matched-haar-information-v1'/ds/'likelihood-diagnostics.npz') as f:haar=f['log_ratios']
        for c in j['cases']:
            records=[b['levels'][-1] for b in c['banks']];d=c['between_bank_logratio_difference'];dif.append(d)
            e=np.array([r['ess'] for r in records]);ess.append(e)
            density_diff.append(np.array(records[0]['log_integrals'])-records[1]['log_integrals'])
            pp=[b['levels'][-1]['log_ratio']-b['levels'][0]['log_ratio'] for b in c['banks']];prefix.append(pp)
            rows.append(dict(dataset=ds,position=c['position'],image_index=c['index'],state=c['state'],
                logratio_bank0=records[0]['log_ratio'],logratio_bank1=records[1]['log_ratio'],bank_difference=d,
                ess_bank0_null=e[0,0],ess_bank0_alternative=e[0,1],ess_bank1_null=e[1,0],ess_bank1_alternative=e[1,1],
                prefix_change_bank0=pp[0],prefix_change_bank1=pp[1],
                haar_32768_logratio_bank0=float(haar[0,-1,1,c['state'],c['index']]),
                haar_32768_logratio_bank1=float(haar[1,-1,1,c['state'],c['index']]),
                maximum_importance_weight=max(max(r['maximum_weight']) for r in records)))
            for i,m in enumerate(c['proposal']['modes']):fitrows.append(dict(dataset=ds,position=c['position'],state=c['state'],mode=i,**m))
        ess=np.array(ess);dif=np.array(dif);density_diff=np.array(density_diff);prefix=np.array(prefix)
        med=np.median(ess,axis=0)
        table.append(f"| {ds} | {np.sum(dif<=.01)}/128 ({np.mean(dif<=.01):.4f}) | {med[0,0]:.1f} / {med[0,1]:.1f} / {med[1,0]:.1f} / {med[1,1]:.1f} | {'pass' if j['gate']['passed'] else 'fail'} | {j['seconds']/60:.2f} |")
        secondary.append(f"| {ds} | {np.max(dif):.6g} | {np.sqrt(np.mean(dif*dif)):.6g} | {np.sqrt(np.mean(density_diff*density_diff)):.6g} | {np.sqrt(np.mean(prefix*prefix)):.6g} | {j['optimizer_successes']}/{j['optimizer_attempts']} |")
        ax=axes[0,col];ax.plot(np.arange(128),dif,'.',ms=4);ax.axhline(.01,color='red',ls='--',lw=1);ax.set(yscale='log',title=f'EMPIAR {ds}',xlabel='Observed image (null/altered paired)',ylabel='Absolute bank log-ratio difference');ax.grid(alpha=.2)
        for bank,style in [(0,'-'),(1,'--')]:
            for model,color in [(0,'#2166ac'),(1,'#b2182b')]:axes[1,col].plot(np.arange(1,129)/128,np.sort(ess[:,bank,model]),ls=style,color=color,label=f'Bank {bank}, model {model}')
        axes[1,col].axhline(256,color='black',lw=1);axes[1,col].set(yscale='log',xlabel='Empirical quantile',ylabel='Importance ESS');axes[1,col].grid(alpha=.2)
    fig.legend(*axes[1,0].get_legend_handles_labels(),loc='outside lower center',ncol=4,fontsize=9)
    for ext in ['png','pdf']:fig.savefig(out/f'integration-gate.{ext}',dpi=180)
    plt.close(fig)
    with (out/'image-diagnostics.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'optimizer-outcomes.json').write_text(json.dumps(fitrows,indent=2)+'\n')
    passes=all(j['gate']['passed'] for _,j,_ in datasets)
    lines=['# Adaptive orientation integration: first frozen gate','','1 October 2026 UTC. All 384 observed-image cases complete: 64 fixed source indices × unchanged/25%-deleted image × three stacks. Every case uses eight optimized proposal modes and two independent 8,192-draw importance banks. The source indices, tolerance and ESS requirements were [declared before outcomes](ADAPTIVE-POSE-INTEGRATION-PROTOCOL.md). All original maps, simulator assumptions and images remain unchanged.','','**Overall gate: '+('PASS.' if passes else 'FAIL.')+'** Passing requires at least 90% of images within .01 log-ratio units between banks and all four median ESS values at least 256 on every stack. No stack is omitted.','','| EMPIAR | Images within tolerance | Median ESS: bank0 null / alt / bank1 null / alt | Stack gate | Minutes |','|---|---:|---:|---|---:|',*table,'','Times are concurrent wall times on the M4 Pro, including image-specific physical Fourier evaluations. They cannot be interpreted as speedups over the cached uniform-template matrix calculation, which shares templates across images. No GPU was rented.','','## Error and optimizer diagnostics','','| EMPIAR | Maximum ratio discrepancy | RMS ratio discrepancy | RMS log-integral discrepancy | RMS last-prefix ratio change | Successful local fits |','|---|---:|---:|---:|---:|---:|',*secondary,'','The same proposal is shared by an image’s two integration banks. Agreement can miss a common uncovered mode; [the numerical requirements note](NUMERICAL-LIKELIHOOD-REQUIREMENTS.md) gives an elementary counterexample. Effective sample size is a weight-concentration diagnostic, not a coverage guarantee. Neither a successful BFGS status nor a positive local Hessian establishes global pose recovery. All 3,072 local-fit outcomes, including nonconvergence, duplicate modes and any reverted terminal point, are retained. No fitted proposal covariance is presented as a pose confidence region.','','## Independent replay','']
    for ds,j,v in datasets:
        lines.append(f"- {ds}: all 128 proposal densities and integral summaries replay. Maximum log-density discrepancy {max(r['density_error'] for r in v['images']):.3g}; maximum log-integral discrepancy {max(r['log_integral_error'] for r in v['images']):.3g}; six direct cell-sum residual checks differ by at most {max(r['error'] for r in v['physical_checks']):.3g}. The prespecified gate verdict reproduces.")
    lines.extend(['','The independent proposal check uses full transformed four-dimensional Gaussian covariance coordinates; it does not call the implemented relative-quaternion density. The physical check explicitly sums all 64³ cells at all retained frequencies for six selected orientations per stack, covering both candidate maps. These checks establish consistency of the computation at those points, not a bound on the remaining orientation integral.','','## Decision boundary','','The first numerical gate does not by itself create an uncertainty method or address measured nonuniform views, image-specific CTFs, amplitude, colored noise, map error, regional specificity or calibration repetitions. The moment branch remains frozen. Under the review-4 stopping rule, at most one substantive integration revision is permitted, directed at the measured failure, with this same cohort, draw count and gate. A second missed prediction ends this numerical branch. A focused independent methodological consultation will inform whether that revision and a specific statistical estimand are worth pursuing. No new full-paper acceptance verdict is claimed.','','All per-image results and importance points are preserved under `results/uncertainty/development/adaptive-pose-integration-v1`; the compact CSV also retains the original uniform-bank scores at the same images.'])
    (ROOT/'research/uncertainty/ADAPTIVE-POSE-INTEGRATION-RESULTS.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main()
