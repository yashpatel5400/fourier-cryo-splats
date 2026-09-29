#!/usr/bin/env python3
"""Derive revision tables/figures only from completed diagnostic records."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    out=BASE/'revision-diagnostics-summary';out.mkdir(exist_ok=True)
    tables=ROOT/'paper/tables';figures=ROOT/'paper/figures';sources={};rows=[];experimental=[]
    tables.mkdir(exist_ok=True);figures.mkdir(exist_ok=True)
    for dataset in ['10028','10049','10076']:
        path=BASE/'experimental-noise-grouped'/f'{dataset}.json';r=json.loads(path.read_text())
        if not r['complete']:raise RuntimeError(f'Incomplete experimental diagnostic: {dataset}')
        sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest();experimental.append(r)
        for record in r['records']:
            rows.append({'dataset':dataset,**{k:v for k,v in record.items() if k!='fit'}})
    with (out/'experimental-cases.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    fig,axes=plt.subplots(1,3,figsize=(7.1,2.55),sharey=True)
    tex=[r'\begin{table*}[t]',r'\centering\small',r'\begin{tabular}{lrrrr}',r'\toprule',
         r'Stack & $n_{\rm inf}/n_{\rm cal}$ & $\sigma_\ell$ (\AA) & Center $\pm$ half-width & Ref.\\',r'\midrule']
    for ax,r in zip(axes,experimental):
        widths=sorted({row['sigma_A'] for row in r['records']},reverse=True)
        for width,color,marker in zip(widths,['#246a9b','#b54c36'],['o','s']):
            data=[row for row in r['records'] if row['sigma_A']==width]
            ax.plot([x['rotation_radius_degrees'] for x in data],[x['relative_half_width'] for x in data],
                    color=color,marker=marker,label=rf'$\sigma_\ell={width:.1f}$ Å')
        ax.set_title('EMPIAR-'+r['dataset'],fontsize=9);ax.set_xlabel('Rotation radius (degrees)',fontsize=8)
        ax.set_xticks([0,1,2]);ax.set_ylim(0,1.05);ax.tick_params(labelsize=8)
        ax.grid(alpha=.18);ax.legend(fontsize=7,frameon=False,loc='lower right')
        broad=[x for x in r['records'] if x['sigma_A']==max(widths) and x['rotation_radius_degrees']==1][0]
        tex.append(f"{r['dataset']} & {r['inference_particles']}/{r['calibration_groups']} & {broad['sigma_A']:.1f} & ${broad['interval_center']:.2f}\\pm{broad['interval_half_width']:.2f}$ & {broad['approximate_reference_target']:.2f}\\\\")
    axes[0].set_ylabel('Half-width / no-data bound',fontsize=8)
    fig.tight_layout();fig.savefig(figures/'revision-experimental.pdf');fig.savefig(out/'experimental-widths.png',dpi=170);plt.close(fig)
    tex.extend([r'\bottomrule',r'\end{tabular}',
       r'\caption{Exploratory experimental sensitivity intervals for the broad central Gaussian average at $1^\circ$ and 0.5\,\AA\ translation radius. Density units use the separately amplitude-aligned unit-norm pilot. One particle per inference/calibration exposure group is used; remaining new groups supply diagnostics. None excludes zero. Ref. is a pilot-pool amplitude-aligned deposited-map value, not a coverage label. Noise commonality, consensus-pose independence and density/pose radii remain unverified.}',
       r'\label{tab:revisionexperimental}',r'\end{table*}'])
    (tables/'revision-experimental.tex').write_text('\n'.join(tex)+'\n')
    path=BASE/'higher-band-ctf-sensitivity.json';ctf=json.loads(path.read_text())
    if not ctf['complete']:raise RuntimeError('Incomplete CTF diagnostic')
    sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    tex=[r'\begin{table*}[t]',r'\centering\small',r'\begin{tabular}{lrr}',r'\toprule',
         r'CTF sensitivity & Relative half-width & No-data cases\\',r'\midrule']
    ctf_summary=[]
    for key,label in [('nominal','Nominal'),('defocus100A',r'Defocus $\pm100$\,\AA'),('defocus500A',r'Defocus $\pm500$\,\AA'),('mixed','Combined nuisance bounds')]:
        selected=[r for r in ctf['records'] if r['setting']==key];width=[r['selected_relative_half_width'] for r in selected]
        fallbacks=sum(r['uses_no_data'] for r in selected)
        tex.append(f'{label} & {min(width):.3f}--{max(width):.3f} & {fallbacks}/6\\\\')
        ctf_summary.append({'setting':key,'min_width':min(width),'max_width':max(width),'no_data':fallbacks})
    tex.extend([r'\bottomrule',r'\end{tabular}',
      r'\caption{Declared CTF-error envelopes added to all six fixed-pose radius-12 fine-feature fits. The combined setting is 100\,\AA\ defocus, $2^\circ$ astigmatism angle, $1^\circ$ phase, 5\% gain and 10\,\AA$^2$ relative B-factor. Radii are sensitivities, not calibrated errors. The triangle bound discards cancellation; its vacuity is not an impossibility result for other methods.}',r'\label{tab:revisionctf}',r'\end{table*}'])
    (tables/'revision-ctf.tex').write_text('\n'.join(tex)+'\n')
    summary={'stage':'Derived revision diagnostics, no new experimental coverage claim','source_sha256':sources,
             'experimental_cases':len(rows),'experimental_excludes_zero':sum(r['excludes_zero'] for r in rows),
             'experimental_reference_inside':sum(r['approximate_reference_inside_interval'] for r in rows),
             'ctf_summary':ctf_summary}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
