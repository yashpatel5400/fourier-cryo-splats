#!/usr/bin/env python3
"""Report every completed independent directional-noise sensitivity case."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    out=BASE/'directional-noise-summary';out.mkdir(exist_ok=True);sources={};studies=[];rows=[]
    for dataset in ['10028','10049','10076']:
        path=BASE/'directional-noise-audit'/f'{dataset}.json';d=json.loads(path.read_text())
        if not d.get('complete'):raise RuntimeError(f'Incomplete source: {dataset}')
        sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest();studies.append(d)
        rows.extend({'dataset':dataset,**{k:v for k,v in row.items() if not isinstance(v,dict)}} for row in d['records'])
    with (out/'cases.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    fig,axes=plt.subplots(2,3,figsize=(7.1,3.8),sharex=True,sharey=True)
    tex=[r'\begin{table*}[t]',r'\centering\small',r'\begin{tabular}{lrrrr}',r'\toprule',
         r'Stack & New cal. $n$ & SD ratio & Center $\pm$ half-width & Ref.\\',r'\midrule']
    for j,d in enumerate(studies):
        widths=sorted({r['sigma_A'] for r in d['records']},reverse=True)
        for i,width in enumerate(widths):
            ax=axes[i,j];selected=[r for r in d['records'] if r['sigma_A']==width];x=[r['rotation_radius_degrees'] for r in selected]
            ax.plot(x,[r['old_relative_half_width'] for r in selected],'o--',color='#989898',label='Original trace')
            ax.plot(x,[r['relative_half_width'] for r in selected],'s-',color='#246a9b',label='Directional')
            ax.set_title(f"{d['dataset']}, $\\sigma_\\ell={width:.1f}$ Å",fontsize=8)
            ax.set_xticks([0,1,2]);ax.set_ylim(0,1.05);ax.grid(alpha=.2);ax.tick_params(labelsize=8)
            if i==1:ax.set_xlabel('Rotation radius (degrees)',fontsize=8)
            if j==0:ax.set_ylabel('Half-width / no-data',fontsize=8)
        broad=next(r for r in d['records'] if r['sigma_A']==widths[0] and r['rotation_radius_degrees']==1)
        tex.append(f"{d['dataset']} & {d['new_calibration_groups']} & {broad['noise_sd_ratio']:.3f} & ${broad['interval_center']:.2f}\\pm{broad['interval_half_width']:.2f}$ & {broad['approximate_reference_target']:.2f}\\\\")
    axes[0,0].legend(fontsize=7,frameon=False);fig.tight_layout()
    fig.savefig(ROOT/'paper/figures/directional-noise.pdf');fig.savefig(out/'widths.png',dpi=170);plt.close(fig)
    tex.extend([r'\bottomrule',r'\end{tabular}',
      r'\caption{Independent directional-noise re-audit of the same broad experimental estimators at $1^\circ$ and 0.5\,\AA. SD ratio compares the new noise bound with the original first-pool trace bound; the new second-pool sample differs. Centers and bias bounds are unchanged before no-data selection. None excludes zero at this pose budget. Ref. remains an approximate deposited-map diagnostic.}',
      r'\label{tab:directionalnoise}',r'\end{table*}'])
    (ROOT/'paper/tables/directional-noise.tex').write_text('\n'.join(tex)+'\n')
    summary={'complete':True,'source_sha256':sources,'cases':len(rows),
      'noise_sd_ratio_range':[min(r['noise_sd_ratio'] for r in rows),max(r['noise_sd_ratio'] for r in rows)],
      'relative_width_range':[min(r['relative_half_width'] for r in rows),max(r['relative_half_width'] for r in rows)],
      'excludes_zero_cases':[{k:r[k] for k in ['dataset','sigma_A','rotation_radius_degrees']} for r in rows if r['excludes_zero']],
      'no_data_cases':sum(r['uses_no_data'] for r in rows),
      'reference_inside_count':sum(r['approximate_reference_inside_interval'] for r in rows)}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
