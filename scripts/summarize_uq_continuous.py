#!/usr/bin/env python3
"""Summarize completed continuous audits and explicit assumption stresses."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development';OUT=BASE/'continuous-summary';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'continuous-audit-v1'})
colors=['#2166ac','#d95f02','#1b9e77'];datasets=['10028','10049','10076']
def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=200,bbox_inches='tight');fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
summary={'stage':'development; same declared continuous cube class across dense and quadrature solvers','comparisons':[]}
fig,axes=plt.subplots(1,2,figsize=(7.2,2.7))
for idx,(dataset,color) in enumerate(zip(datasets,colors)):
    dense=json.loads((BASE/'continuous-optimized'/f'{dataset}.json').read_text());quad=json.loads((BASE/'continuous-quadrature-optimized'/f'{dataset}.json').read_text())
    assert len(dense['targets'])==len(quad['targets'])==4
    values=[r['width_fraction_of_no_data'] for r in quad['targets']];axes[0].bar(np.arange(4)+(idx-1)*.23,values,.22,label=dataset,color=color)
    cell=json.loads((BASE/'continuous-cube_cells'/f'{dataset}.json').read_text());record=next(r for r in cell['targets'] if r['target']=='center' and r['width_fraction_field']==.07)
    axes[1].plot([r['box'] for r in record['cell_subspaces']],[r['continuous_boundary_coverage_using_cell_width'] for r in record['cell_subspaces']],'-o',color=color,markersize=3)
    for d,q in zip(dense['targets'],quad['targets']):
        assert (d['target'],d['width_fraction_field'])==(q['target'],q['width_fraction_field'])
        summary['comparisons'].append({'dataset':dataset,'target':d['target'],'width':d['width_fraction_field'],
            'relative_width':q['width_fraction_of_no_data'],'dense_quadrature_relative_difference':abs(q['fit']['half_width']/d['fit']['half_width']-1),
            'continuous_gap':q['fit']['history'][-1]['relative_gap'],'reference_coverage':q['reference_analytic_coverage']})
axes[0].set_xticks(range(4),['Center\n.03','Contrast\n.03','Center\n.07','Contrast\n.07']);axes[0].set_ylabel('Continuous width / no-data width');axes[0].set_xlabel('Target and width / field');axes[0].set_ylim(0,1)
axes[1].set_xlabel('Cell grid used for the restricted audit');axes[1].set_ylabel('Coverage at continuous\nbias boundary');axes[1].set_ylim(0,1.03);axes[1].axhline(.95,color='black',ls=':',lw=1)
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,ncol=3,frameon=False,loc='upper center',title='EMPIAR acquisition geometry');fig.tight_layout(rect=[0,0,1,.82]);save(fig,'continuous-audit')

fig,axes=plt.subplots(2,3,figsize=(7.2,4.5),sharey=True)
settings=[('noise_scale','boundary_for_fixed_pose','Gaussian noise scale / assumed'),('noise_correlation','boundary_for_fixed_pose','Within-particle noise correlation'),
          ('support_violation','outside_direction_for_fixed_pose','Outside-support energy fraction'),('two_state_heterogeneity','view_estimator_coupled_for_fixed_pose','Two-state density amplitude'),
          ('joint_pose_radius','shared_estimator_nonlinear_adversary','Actual / assumed pose radius'),('coherent_defocus_error_A','independent_map','Coherent defocus error (A)')]
for dataset,color in zip(datasets,colors):
    source=json.loads((BASE/'assumption-stress'/f'{dataset}-center.json').read_text())
    for ax,(axis,signal,xlabel) in zip(axes.ravel(),settings):
        for method,style in [('fixed_pose','--'),('quadratic_shared','-')]:
            rows=sorted([r for r in source['records'] if r['axis']==axis and r['signal']==signal and r['method']==method],key=lambda r:r['setting'])
            ax.plot([r['setting'] for r in rows],[r['analytic_coverage'] for r in rows],style,color=color,marker='o',markersize=2,label=dataset if style=='-' else None)
        ax.set_xlabel(xlabel);ax.axhline(.95,color='black',ls=':',lw=.6);ax.set_ylim(-.03,1.04)
axes[0,0].set_ylabel('Analytic coverage');axes[1,0].set_ylabel('Analytic coverage');handles,labels=axes[0,0].get_legend_handles_labels()
fig.legend(handles,labels,ncol=3,frameon=False,loc='upper center',title='Dashed: fixed pose; solid: finite-voxel shared-density pose audit')
fig.tight_layout(rect=[0,0,1,.89]);save(fig,'assumption-stress-center')
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(OUT/'README.md').write_text('''# Continuous and failure-control development figures

The continuous density class is the full unit cube, whereas the earlier pose
experiments use a supported finite voxel space. Do not compare their widths as
if the support assumption were identical. Continuous bias-boundary coverage
uses an estimator-specific allowed adversary, not unknown experimental truth.

The failure-control plot uses central targets and a common fixed-pose boundary
signal within each geometry for noise tests. Support and heterogeneous controls
use the fixed-pose estimator's selected direction; pose controls use the saved
nonlinear pose adversary. Each control is shared between displayed estimators.
The global estimator-aligned artificial noise mode and other signals/settings
are retained in the full result JSONs. These are deliberate assumption stresses,
not representative estimates of failure rates in experimental data.
''')
