#!/usr/bin/env python3
"""Summarize completed development coverage without pseudo-independent trials."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'results/uncertainty/development/coverage-benchmark.json'
out=source.parent/'coverage-summary';out.mkdir(exist_ok=True)
data=json.loads(source.read_text());records=[]
for case in data['cases']:
    for row in case['rows']:
        if row['truth_case']=='bounds_violated':continue
        records.append({**{k:case[k] for k in ['seed','geometry','angle_deg']},**row})
frame=pd.DataFrame(records);summary=[]
for (geometry,angle,method),part in frame.groupby(['geometry','angle_deg','method']):
    # Width is constant for a fixed design/method; count each design once.
    per_design=part.drop_duplicates('seed')
    summary.append({'geometry':geometry,'angle_deg':angle,'method':method,'independent_designs':len(per_design),
                    'signal_design_pairs':len(part),'minimum_analytic_coverage':part.analytic_coverage.min(),
                    'mean_analytic_coverage':part.analytic_coverage.mean(),
                    'pooled_mc_coverage_descriptive_only':part.covered.sum()/part.noise_replicates.sum(),
                    'median_width_fraction_of_no_data':per_design.width_fraction_of_no_data.median(),
                    'minimum_width_fraction_of_no_data':per_design.width_fraction_of_no_data.min(),
                    'maximum_width_fraction_of_no_data':per_design.width_fraction_of_no_data.max()})
table=pd.DataFrame(summary);table.to_csv(out/'stratified.csv',index=False,lineterminator='\n')
methods=['posterior','joint_gaussian','fixed_pose_bias_aware','linearized_nuisance','coarse_remainder','structured_interaction']
labels=['Gaussian posterior','Joint Gaussian pose','Fixed-pose bias-aware','Linearized pose','Coarse remainder','Structured interaction']
colors=['#c54b43','#a4659b','#d99736','#668cb8','#777777','#167c75']
plt.rcParams.update({'font.size':8,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':7,'svg.hashsalt':'fourier-cryo-splats-uq'})
fig,axes=plt.subplots(2,2,figsize=(7,4.8),layout='constrained')
for col,geometry in enumerate(['uniform','preferred']):
    for i,(method,label,color) in enumerate(zip(methods,labels,colors)):
        part=table[(table.geometry==geometry)&(table.method==method)].sort_values('angle_deg')
        x=np.arange(len(part))+(i-2.5)*.11
        axes[0,col].scatter(x,part.minimum_analytic_coverage,s=20,label=label,color=color)
        axes[1,col].scatter(x,part.median_width_fraction_of_no_data,s=20,color=color)
    for row in range(2):
        ax=axes[row,col];ax.set_xticks([0,1],['0.5°','2°']);ax.set_xlabel('Rotation uncertainty radius');ax.grid(axis='y',alpha=.2)
    axes[0,col].axhline(.95,color='black',ls='--',lw=1)
    axes[1,col].axhline(1,color='black',ls='--',lw=1)
    axes[0,col].set_ylim(-.04,1.04);axes[1,col].set_ylim(0,1.1);axes[0,col].set_title(geometry.capitalize()+' views')
axes[0,0].set_ylabel('Minimum analytic coverage\nover tested bounded signals')
axes[1,0].set_ylabel('Median interval width / no-data width')
fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False)
fig.suptitle('Development: exact finite dictionary, known bounds, fixed pilot',fontsize=9)
fig.savefig(out/'coverage-width.png',dpi=240);fig.savefig(out/'coverage-width.svg',metadata={'Date':None});plt.close(fig)
svg=out/'coverage-width.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
notes='''# Development coverage summary

Source: `../coverage-benchmark.json`; 48 independent design settings, with paired
noise draws across methods within each signal. These are exact-dictionary,
fixed-pilot, known-bound simulations. They are not real-data calibration.

Coverage is conditional on each fixed signal; minimum means minimum over the
finite tested set, not a proof of the worst case. Random bounded and separately
constructed adversarial signals are both included. The bounds-violated stress
is excluded from this table. Widths count each independent design once.

The pooled Monte Carlo proportion is descriptive. Noise draws are paired across
methods, and multiple signals share a design, so treating all rows as independent
binomial trials would give invalid comparison error bars. Full Gaussian
posteriors separately passed the exactly matched prior-predictive check at 0.95;
fixed-parameter boundary failures do not contradict that result.

The unregularized baseline is present in the CSV (including its very large
widths); the figure focuses on the six main methods to remain legible. The
incomplete run is not included. This is an exploratory, not confirmatory, study.
'''
(out/'README.md').write_text(notes)
print(table[table.method.isin(['structured_interaction','linearized_nuisance'])].to_string(index=False))
