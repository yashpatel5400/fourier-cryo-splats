#!/usr/bin/env python3
"""Generate the post-review-3 tables from complete saved records."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'results/uncertainty/development/refitting-bias-summary-v1/summary.json').read_text())
lookup={tuple(r[k] for k in ['dataset','template','target','method','image_mode']):r for r in d['groups']}
lines=[r'\begin{table*}[t]',r'\centering\footnotesize',r'\caption{Post hoc audit-weight diagnostics. Each pair is same-image / independent-image; every cell has 200 paired trials. RMSE is divided by the fixed pilot absolute error. The complete release includes all five distinct estimators, all radius choices, and individual Monte Carlo intervals.}',r'\label{tab:refittingbias}',r'\begin{tabular}{lllrrr}\toprule',r'Stack & Pose control & Target & Mean error / SD & Noise-only coverage & RMSE / pilot error \\\midrule']
for ds in ['10028','10049','10076']:
 for template,label in [('true_pose','True poses'),('oracle_reference','Oracle template'),('independent_pilot','Pilot template')]:
  for target in ['center','contrast']:
   a=lookup[ds,template,target,'fixed_folded','same_image'];b=lookup[ds,template,target,'fixed_folded','independent_image']
   lines.append(f"{ds} & {label} & {target} & {a['standardized_error_mean']:.3f} / {b['standardized_error_mean']:.3f} & {a['B0_coverage']:.3f} / {b['B0_coverage']:.3f} & {a['rmse_over_pilot_error']:.3f} / {b['rmse_over_pilot_error']:.3f} "+r'\\')
lines += [r'\bottomrule\end{tabular}',r'\end{table*}']
(ROOT/'paper/tables/refitting-bias-detail.tex').write_text('\n'.join(lines)+'\n')
rows=[]
for ds in ['10028','10049','10076']:
 records=json.loads((ROOT/f'results/uncertainty/development/continuous-gaussian-review2-v2/{ds}.json').read_text())['records']
 for tau in [2.,1.]:
  for mode in ['fixed_pose','local_gaussian_pose']:
   values=[r['fit']['fixed_pose_uniform_class_coverage'] for r in records if r['prior_directional_sd']==tau and r['pose_model']==mode]
   assert len(values)==4
   rows.append(dict(dataset=ds,prior_scale=tau,pose_model=mode,coverage_min=min(values),coverage_max=max(values),fits=4))
s=r'''\begin{table}[t]
\centering\small
\caption{Minimum fixed-pose uniform class coverage across four targets, at nominal $1-.05/12=.995833$. Columns separate the prior scale and fixed/linearized covariance weights. These are analytic class minima, not coverage on a selected generator.}
\label{tab:gpclass}
\begin{tabular}{lrrrr}\toprule
Stack & $2$, fixed & $2$, linear & $1$, fixed & $1$, linear \\\midrule
'''
for ds in ['10028','10049','10076']:
 vals=[r['coverage_min'] for r in rows if r['dataset']==ds];s+=ds+' & '+' & '.join(f'{v:.6f}' for v in vals)+r' \\'+'\n'
s+=r'\bottomrule\end{tabular}\end{table}'+'\n'
(ROOT/'paper/tables/continuous-gaussian-class-coverage.tex').write_text(s)
(ROOT/'results/uncertainty/development/refitting-bias-summary-v1/gaussian-class-coverage.json').write_text(json.dumps(rows,indent=2)+'\n')
