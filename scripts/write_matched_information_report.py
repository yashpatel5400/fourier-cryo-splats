#!/usr/bin/env python3
"""Report every frozen stack and quadrature level without selecting a new test."""
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
    out=BASE/'matched-information-report-v2'
    if out.exists():raise ValueError('Preserve report attempts')
    out.mkdir()
    check=json.loads((BASE/'matched-information-verification-v1/summary.json').read_text())
    assert check['complete']
    text=['# Matched information and nuisance diagnosis', '',
        '1 October 2026 UTC. All three frozen stage-A/stage-B comparisons complete. '
        'These are post hoc diagnostics on the same candidate, 20 Å region, 220 complex frequencies, '
        'transfer profile and 8,192 paired images per stack. The [protocol](MATCHED-INFORMATION-LEDGER-PROTOCOL.md) '
        'preceded their outcomes. The original simulator still assumes Haar orientations, one CTF, '
        'amplitude one and known independent unit real/imaginary noise. No experimental calibration or '
        'new test is claimed.', '',
        '## Matched separation at 25% deletion', '',
        'Define D(T) = (E_alt T − E_null T)² / Var_null T. Moment means/conditional variances use '
        'exact Gaussian formulas, with Monte Carlo integration over views. The known-pose column is '
        'twice Gaussian KL averaged over these views; it also equals D for a conditionally centered '
        'known-pose linear score. The numerical likelihood columns are D of approximate Haar mixture '
        'log ratios, not certified KL values or information ceilings. Sampling standard errors for the '
        'paired likelihood mean gap and known-pose energy are retained in the full JSON; '
        'they do not include quadrature bias.', '',
        '| EMPIAR | Known-pose 2 KL | Power matched / Fisher D | Combined matched / Fisher D | Numerical likelihood D, banks 0 / 1 |',
        '|---|---:|---:|---:|---:|']
    rows=[]; stability=[]; nuisance=[];components=[]
    fig,axes=plt.subplots(1,3,figsize=(12,3.6),constrained_layout=True)
    for ax,ds in zip(axes,['10028','10049','10076']):
        a=json.loads((BASE/'matched-information-ledger-v1'/ds/'summary.json').read_text())
        b=json.loads((BASE/'matched-haar-information-v1'/ds/'summary.json').read_text())
        assert a['complete'] and b['complete']
        cases={c['key']:c for c in a['cases']}
        def dc(key):return next(c for c in cases[key]['deletions'] if c['deletion']==.25)
        known=dc('power_matched')['known_pose_twice_kl']
        selected=[c for c in b['cases'] if c['deletion']==.25 and c['quadrature_views']==32768]
        text.append(f"| {ds} | {known:.6g} | {dc('power_matched')['squared_separation']:.6g} / {dc('power_fisher')['squared_separation']:.6g} | {dc('power_bispectrum_matched')['squared_separation']:.6g} / {dc('power_bispectrum_fisher')['squared_separation']:.6g} | {selected[0]['squared_separation']:.6g} / {selected[1]['squared_separation']:.6g} |")
        with np.load(BASE/'matched-haar-information-v1'/ds/'likelihood-diagnostics.npz') as f:
            lr=f['log_ratios']
        for di,delta in enumerate(b['deletions']):
            for state in [0,1]:
                x,y=lr[0,-1,di,state],lr[1,-1,di,state]
                stability.append(dict(dataset=ds,deletion=delta,state=state,
                    bank_rms=float(np.sqrt(np.mean((x-y)**2))),bank_mean_difference=float(np.mean(x-y)),
                    bank_correlation=float(np.corrcoef(x,y)[0,1]),score_sd_bank0=float(x.std(ddof=1)),
                    bank0_last_doubling_rms=float(np.sqrt(np.mean((x-lr[0,-2,di,state])**2))),
                    bank1_last_doubling_rms=float(np.sqrt(np.mean((y-lr[1,-2,di,state])**2)))))
        for c in b['cases']:
            d=c['diagnostics'][0]
            rows.append(dict(dataset=ds,**{k:c[k] for k in ['bank','quadrature_views','deletion','null_logratio_mean','alternative_logratio_mean','signed_mean_gap','paired_gap_standard_error','null_variance','squared_separation','known_pose_twice_kl']},
                null_ess_median=d['ess_quantiles']['0.5'],null_max_weight_q95=d['maximum_weight_quantiles']['0.95']))
        for key,c in cases.items():
            z=dc(key)
            if key.endswith('_component'):
                components.append(f"| {ds} | {key} | {z['squared_separation']:.6g} |")
            else:
                nuisance.append(f"| {ds} | {key} | {c['power_weight_sum']:.6g} | {100*z['white_noise_variance_error_matching_gap']:.4g}% | {100*c['view_mean_variance']/c['null_total_variance']:.3g}% | {'yes' if z['outside_global_amplitude_mean_envelope'] else 'no'} / {'yes' if z['outside_view_dependent_mean_envelope'] else 'no'} |")
        for bank in [0,1]:
            cc=[c for c in b['cases'] if c['bank']==bank and c['deletion']==.25]
            ax.plot([c['quadrature_views'] for c in cc],[c['squared_separation'] for c in cc],marker='o',label=f'Mixture bank {bank}')
        for key,label,style,color in [('power_matched','Power matched',':','#1a9850'),('power_fisher','Power Fisher','--','#1a9850'),('power_bispectrum_matched','Combined matched',':','#762a83'),('power_bispectrum_fisher','Combined Fisher','--','#762a83')]:
            ax.axhline(dc(key)['squared_separation'],ls=style,lw=1,label=label,color=color)
        ax.axhline(known,color='black',lw=1.1,label='Known-pose 2 KL')
        ax.set(xscale='log',yscale='log',title=f'EMPIAR {ds}',xlabel='Quadrature views per bank',ylabel='Squared mean separation D')
        ax.set_xticks(b['levels'],['4k','8k','16k','32k'])
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator());ax.grid(alpha=.2)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=4,fontsize=8)
    for ext in ['png','pdf']:fig.savefig(out/f'matched-separation.{ext}',dpi=180)
    plt.close(fig)
    for name,data in [('quadrature-comparison',rows),('quadrature-stability',stability)]:
        with (out/(name+'.csv')).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    text.extend(['','## Integration is not converged','','| EMPIAR | Null ESS median, bank 0 / 1 | Bank RMS / score SD, null | Last doubling RMS, bank 0 / 1 |','|---|---:|---:|---:|'])
    for ds in ['10028','10049','10076']:
        cc=[r for r in rows if r['dataset']==ds and r['deletion']==.25 and r['quadrature_views']==32768]
        s=next(r for r in stability if r['dataset']==ds and r['deletion']==.25 and r['state']==0)
        text.append(f"| {ds} | {cc[0]['null_ess_median']:.4g} / {cc[1]['null_ess_median']:.4g} | {s['bank_rms']:.5g} / {s['score_sd_bank0']:.5g} | {s['bank0_last_doubling_rms']:.5g} / {s['bank1_last_doubling_rms']:.5g} |")
    text.extend(['', 'The mean separation agrees reasonably between banks, while individual log ratios and the '
        'last quadrature doubling remain materially unstable. In particular, typical 10028 images have '
        'nearly single-orientation support. A shared ratio can partly cancel density errors, but that '
        'does not certify accurate mixture densities or calibrated inference. More uniform random '
        'orientations alone is not the next statistical-method variant; the unresolved numerical '
        'integration must be addressed directly. These diagnostics justify investigating a likelihood '
        'approach, not reporting its asymptotic efficiency or treating failure as impossibility.', '',
        '## Amplitude and noise sensitivity', '',
        '| EMPIAR | Frozen score | Sum of power weights | White-noise variance error matching deletion gap | View share of null variance | Alternative outside amplitude envelope: global / view-dependent |',
        '|---|---|---:|---:|---:|---|',*nuisance,'',
        'The amplitude interval is [.9,1.1]. These are **mean-score** envelopes, not the original '
        'event-probability null bounds. For 10028 and 10076, view-dependent amplitude can cover the '
        '25%-deletion mean gap for every original direction. Even without that ambiguity, roughly '
        '0.3% white-noise variance error matches the 10076 deletion gap. The bispectrum mean has '
        'zero noise shift only under the specified independent nonredundant Gaussian coordinates. '
        'The calculation does not quantify arbitrary colored or correlated noise. The experimental '
        '[metadata inventory](STACK-NUISANCE-INVENTORY-RESULTS.md) does not establish this precision.', '',
        '## Existing combined score components', '',
        'These are the unchanged coefficients of each combined score, with one block zeroed for '
        'diagnosis. Neither block is refit or promoted to a new test. Their separations do not add '
        'because their noise and view covariances are nonzero.', '',
        '| EMPIAR | Component | D at 25% deletion |','|---|---|---:|',*components,'',
        '## Verification and remaining work','',
        'All 60 original amplitude-one score prefixes replay within 1.01e-11. Independent complex '
        'Gaussian monomial enumeration checks 216 conditional variances across three noise levels '
        'and three fixed image indices; the largest difference is 1.78e-15. Direct Euclidean Gaussian '
        'residuals plus SciPy log-sum-exp reproduce 576 saved likelihood ratios within 8.53e-14, '
        'including quadrature ESS and maximum weights. All saved summary separations are replayed. '
        'These are floating-point checks, not rigorous integration error bounds.', '',
        'The full ledger still lacks the original **event** envelope decomposition, finite-calibration '
        'slack at matched draws, measured experimental noise spectra and realistic-nuisance likelihood '
        'validation. The surviving scientific question is whether a tractable nuisance-aware likelihood '
        'can support a clearly specified uncertainty target. No new score fitting, calibration repetition, '
        'acceptance-level claim or full-review request follows from these numerical comparisons.', '',
        'Sources: `scripts/analyze_matched_information.py`, `scripts/analyze_matched_haar_likelihood.py`, '
        '`scripts/verify_matched_information.py`; complete summaries and raw arrays are under '
        '`results/uncertainty/development/matched-*-v1`. The reporting CSV includes all four deletion '
        'fractions, both banks and every declared prefix. Earlier arrays and failed method branches remain unchanged.',''])
    (ROOT/'research/uncertainty/MATCHED-INFORMATION-LEDGER-RESULTS.md').write_text('\n'.join(text))


if __name__=='__main__':main()
