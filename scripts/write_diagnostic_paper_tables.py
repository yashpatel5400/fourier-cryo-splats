#!/usr/bin/env python3
"""Small manuscript tables generated from verified experiment records."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def write(name,caption,label,columns,header,rows):
    text=['\\begin{table}[t]','\\centering\\small',f'\\caption{{{caption}}}',f'\\label{{{label}}}',
          f'\\begin{{tabular}}{{{columns}}}','\\toprule',header+'\\\\','\\midrule',
          *[r+'\\\\' for r in rows],'\\bottomrule','\\end{tabular}','\\end{table}']
    (ROOT/'paper/tables'/name).write_text('\n'.join(text)+'\n')


def main():
    separation=[];integration=[];population=[]
    for ds in ['10028','10049','10076']:
        moment=json.loads((BASE/'matched-information-ledger-v1'/ds/'summary.json').read_text())
        haar=json.loads((BASE/'matched-haar-information-v1'/ds/'summary.json').read_text())
        keys={'power_matched','power_fisher','power_bispectrum_matched','power_bispectrum_fisher'}
        best=max(d['squared_separation'] for c in moment['cases'] if c['key'] in keys for d in c['deletions'] if d['deletion']==.25)
        h=next(c for c in haar['cases'] if c['bank']==0 and c['quadrature_views']==32768 and c['deletion']==.25)
        separation.append(f"{ds} & {h['known_pose_twice_kl']:.5f} & {best:.6f} & {h['squared_separation']:.5f}")
        j=json.loads((BASE/'catalog-pose-integration-v1'/ds/'summary.json').read_text())
        old=json.loads((BASE/'adaptive-pose-integration-v1'/ds/'summary.json').read_text())
        assert j['complete'] and json.loads((BASE/'catalog-pose-integration-v1'/ds/'independent-check.json').read_text())['complete']
        oldn=sum(c['between_bank_logratio_difference']<=.01 for c in old['cases']);newn=sum(c['between_bank_logratio_difference']<=.01 for c in j['cases'])
        rms=(sum(c['between_bank_logratio_difference']**2 for c in j['cases'])/128)**.5
        integration.append(f"{ds} & {oldn}/128 & {newn}/128 & {rms:.6f}")
        pop=json.loads((BASE/'population-materiality-v1'/ds/'summary.json').read_text())
        assert pop['complete'] and json.loads((BASE/'population-materiality-v1'/ds/'independent-check.json').read_text())['complete']
        pairs=next(g for g in pop['grids'] if g['bins_per_axis']==8)['pairs']
        bias=max(abs(p['bias']) for p in pairs);lo=min(p['amplitude_tangent_harmonic'] for p in pairs);hi=max(p['amplitude_tangent_harmonic'] for p in pairs)
        population.append(f"{ds} & {bias:.6f} & {lo:.3f}--{hi:.3f} & Fail")
    write('diagnostic-separation.tex','Matched 25\\% deletion. Known-pose twice-KL is an oracle quantity; the other columns are squared separation $D$ of the best tested moment and bank-0 template statistics.','tab:diag-separation','lrrr','EMPIAR & Oracle $2\\,\\mathrm{KL}$ & Moment $D$ & Bank $D$',separation)
    write('diagnostic-integration.tex','Unchanged .01 ratio-discrepancy criterion. Every repaired stack also passes the median-ESS requirement.','tab:diag-integration','lrrr','EMPIAR & First & Repair & Repair RMS',integration)
    write('diagnostic-population.tex','Primary population screen: largest absolute bias and harmonic amplitude-tangent surrogate range across all prelisted assignments. The same pair must meet both criteria; all fail materiality.','tab:diag-population','lrrl','EMPIAR & Max. $|t_\\star-p|$ & Harmonic range & Gate',population)


if __name__=='__main__':main()
