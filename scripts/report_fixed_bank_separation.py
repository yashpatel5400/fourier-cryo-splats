#!/usr/bin/env python3
"""Sampling errors of a fixed statistic; no KL, ceiling or power guarantee."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    out=BASE/'fixed-bank-separation-v1'
    if out.exists():raise ValueError('Preserve all report attempts')
    out.mkdir();rows=[];inputs={};table=[]
    for ds in ['10028','10049','10076']:
        directory=BASE/'matched-haar-information-v1'/ds
        j=json.loads((directory/'summary.json').read_text());assert j['complete']
        ap=directory/'likelihood-diagnostics.npz';assert sha(ap)==j['arrays']['sha256']
        inputs[ds]={str(p.relative_to(ROOT)):sha(p) for p in [directory/'summary.json',ap]}
        with np.load(ap) as f:ratios=f['log_ratios']
        for c in j['cases']:
            b=c['bank'];li=j['levels'].index(c['quadrature_views']);di=j['deletions'].index(c['deletion'])
            x,y=ratios[b,li,di];n=len(x);h=y-x;mu=h.mean();variance=x.var(ddof=1);factor=n/(n-1)
            value=mu*mu/variance
            influence=2*mu/variance*(h-mu)-mu*mu/variance**2*(factor*(x-x.mean())**2-variance)
            se=float(influence.std(ddof=1)/np.sqrt(n))
            np.testing.assert_allclose(value,c['squared_separation'],rtol=1e-12,atol=1e-12)
            np.testing.assert_allclose(h.std(ddof=1)/np.sqrt(n),c['paired_gap_standard_error'],rtol=1e-12,atol=1e-12)
            # Independent numerical directional derivative for three empirical
            # contamination directions, with the same fixed finite-n factor.
            gradient_errors=[]
            eps=1e-6
            for index in [0,n//2,n-1]:
                def contaminated(t):
                    mh=(1-t)*mu+t*h[index]
                    mx=(1-t)*x.mean()+t*x[index]
                    second=(1-t)*np.mean(x*x)+t*x[index]**2
                    return mh*mh/(factor*(second-mx*mx))
                derivative=(contaminated(eps)-contaminated(-eps))/(2*eps)
                error=abs(derivative-influence[index]);gradient_errors.append(float(error))
                if error>1e-6*max(1,abs(influence[index])):raise ArithmeticError('Delta influence derivative disagrees')
            row=dict(dataset=ds,bank=b,views=c['quadrature_views'],deletion=c['deletion'],
                     signed_mean_gap=float(mu),paired_gap_se=c['paired_gap_standard_error'],
                     null_variance=float(variance),squared_separation=float(value),
                     delta_sampling_se=se,maximum_influence_derivative_error=max(gradient_errors))
            rows.append(row)
            if row['views']==32768 and row['deletion']==.25:
                table.append(f"| {ds} | {b} | {mu:.6g} ± {row['paired_gap_se']:.3g} | {value:.6g} ± {se:.3g} |")
    with (out/'all-statistics.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'summary.json').write_text(json.dumps(dict(complete=True,cases=len(rows),inputs=inputs,
        source_sha256=sha(Path(__file__)),maximum_influence_derivative_error=max(r['maximum_influence_derivative_error'] for r in rows)),indent=2)+'\n')
    lines=['# Fixed-bank statistic separation and sampling error','',
        '1 October 2026 UTC. Post hoc reporting addition to the existing matched ledger. It uses every saved bank/level/deletion (96 rows), with no new simulated images, orientation draws or score tuning. A finite template-bank log ratio is a statistic even when it is an inaccurate Haar marginal likelihood. Its empirical squared separation D=(E[T1-T0])²/Var(T0) can be described without calling it a KL divergence or information ceiling.',
        '', '| EMPIAR | Bank | Paired mean gap ± MC SE | Squared separation D ± delta MC SE |',
        '|---|---|---:|---:|',*table,'',
        'The table shows the pre-existing 25%-deletion, 32,768-template entries; all 96 entries remain in the CSV. SEs concern 8,192 independent simulated source images conditional on the fixed bank, maps, CTF and unit-white-noise/Haar model. The null and alternative of a source image share noise and pose. Bank errors are correlated and are not independent replications of a clinical or biological finding. No experimental model-error uncertainty is included.',
        '', 'For h_i=T1_i-T0_i, mean gap mu and sample null variance v, the first-order influence used for D is 2 mu/v (h_i-mu) - mu²/v² [n/(n-1)(T0_i-mean(T0))²-v]. Its sample standard deviation divided by sqrt(n) is the reported asymptotic delta SE. This accounts for estimated variance and covariance with the gap. It is not a finite-sample confidence bound or a normal-theory power guarantee. Three independent finite-difference contamination derivatives per row verify the implemented influence. The gap SE alone is not relabelled as an SE of D.',
        '', 'This reporting correction accepts D1 of the [focused consultation](reviews/post-round04-method-consultation/critique.md) only as a descriptive separation comparison. The existing template statistic shows substantially larger separation than the moments under the matched simulator. Its instability still prevents treating it as a converged marginal likelihood, a statistical impossibility bound or a likelihood certified uniformly over parameters. The improved continuous integration gate uses different image-specific proposals and does not retroactively change these frozen template-bank scores.']
    (ROOT/'research/uncertainty/FIXED-BANK-SEPARATION-RESULTS.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main()
