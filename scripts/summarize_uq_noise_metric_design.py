#!/usr/bin/env python3
"""Compare all covariance-guided fits with the same-split directional audit."""
import csv
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    out=BASE/'noise-metric-summary';out.mkdir(exist_ok=True);sources={};rows=[]
    tex=[r'\begin{table*}[t]',r'\centering\small',r'\begin{tabular}{llrrrr}',r'\toprule',
      r'Stack & $\sigma_\ell$ (\AA) & Fixed / old & Fixed / new & $1^\circ$ / old & $1^\circ$ / new\\',r'\midrule']
    for dataset in ['10028','10049','10076']:
        fits={}
        for label,folder in [('old','directional-noise-audit'),('new','noise-metric-design')]:
            path=BASE/folder/f'{dataset}.json';d=json.loads(path.read_text())
            if not d['complete']:raise RuntimeError(f'Incomplete {folder}/{dataset}')
            sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest();fits[label]=d
        lookup={(r['sigma_A'],r['rotation_radius_degrees']):r for r in fits['old']['records']}
        for r in fits['new']['records']:
            old=lookup[r['sigma_A'],r['rotation_radius_degrees']]
            rows.append({'dataset':dataset,'sigma_A':r['sigma_A'],'degrees':r['rotation_radius_degrees'],
              'old_relative_width':old['relative_half_width'],'new_relative_width':r['relative_half_width'],
              'new_uses_no_data':r['uses_no_data'],'new_excludes_zero':r['excludes_zero'],
              'new_reference_inside':r['approximate_reference_inside_interval'],
              'proxy_fit_gap':r['design_proxy_fit']['history'][-1]['relative_gap'],
              'proxy_fit_converged':r['design_proxy_fit']['converged']})
        for width in sorted({r['sigma_A'] for r in fits['new']['records']},reverse=True):
            zero=next(r for r in rows if r['dataset']==dataset and r['sigma_A']==width and r['degrees']==0)
            one=next(r for r in rows if r['dataset']==dataset and r['sigma_A']==width and r['degrees']==1)
            tex.append(f"{dataset} & {width:.1f} & {zero['old_relative_width']:.3f} & {zero['new_relative_width']:.3f} & {one['old_relative_width']:.3f} & {one['new_relative_width']:.3f}\\\\")
    tex.extend([r'\bottomrule',r'\end{tabular}',
      r'\caption{Relative half-widths before and after covariance-guided fixed-pose refitting. Both variants independently calibrate noise on the same second exposure half. A training second-moment proxy guides the new weights; it is not a noise confidence bound. All six fixed-pose widths improve, but five one-degree widths worsen and the sixth remains vacuous. At two degrees every new fit falls back to no data.}',
      r'\label{tab:noisemetric}',r'\end{table*}'])
    (ROOT/'paper/tables/noise-metric-design.tex').write_text('\n'.join(tex)+'\n')
    with (out/'cases.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    summary={'complete':True,'source_sha256':sources,'cases':len(rows),'no_data_cases':sum(r['new_uses_no_data'] for r in rows),
      'excludes_zero_cases':[{k:r[k] for k in ['dataset','sigma_A','degrees']} for r in rows if r['new_excludes_zero']],
      'reference_inside_count':sum(r['new_reference_inside'] for r in rows),
      'largest_proxy_optimization_gap':max(r['proxy_fit_gap'] for r in rows),
      'fixed_pose_improvements':sum(r['new_relative_width']<r['old_relative_width'] for r in rows if r['degrees']==0),
      'one_degree_worsenings':sum(r['new_relative_width']>r['old_relative_width'] for r in rows if r['degrees']==1)}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
