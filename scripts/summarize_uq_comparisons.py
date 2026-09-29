#!/usr/bin/env python3
"""Plot completed development comparisons without selecting favorable settings."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, LogLocator
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development';OUT=BASE/'comparison-summary';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'uq-comparison-20260929'})
DATASETS=['10028','10049','10076'];COLORS=['#2166ac','#d95f02','#1b9e77']


def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=180,bbox_inches='tight')
    path=OUT/(name+'.svg');fig.savefig(path,bbox_inches='tight',metadata={'Date':None})
    path.write_text('\n'.join(s.rstrip() for s in path.read_text().splitlines())+'\n');plt.close(fig)


summary={'stage':'development, all completed settings retained','quadratic_pose':[],'bootstrap':[],'neural':[]}
fig,axes=plt.subplots(2,2,figsize=(7.2,4.8),sharex='col')
for col,target in enumerate(['contrast','center']):
    source=json.loads((BASE/f'ambient-quadratic-{target}.json').read_text());assert len(source['cases'])==3
    for case,color in zip(source['cases'],COLORS):
        angles=[a['angle_deg'] for a in case['angles']]
        for method,style in [('bounded_pose_ambient','--'),('quadratic_pose_ambient','-')]:
            widths=[a['fits'][method]['width_fraction_of_no_data'] for a in case['angles']]
            powers=[next(r['correct_sign_probability'] for r in a['coverage'] if r['method']==method and r['truth']=='independent_EMDB_reference_coherent_pose') for a in case['angles']]
            axes[0,col].plot(angles,widths,style,color=color,marker='o',markersize=3,label=case['dataset'] if style=='-' else None)
            axes[1,col].plot(angles,powers,style,color=color,marker='o',markersize=3)
            for angle,width,power in zip(case['angles'],widths,powers):
                summary['quadratic_pose'].append({'dataset':case['dataset'],'target':target,'angle_deg':angle['angle_deg'],'method':method,
                    'width_fraction':width,'reference_sign_power':power,
                    'minimum_tested_coverage':min(r['analytic_coverage'] for r in angle['coverage'] if r['method']==method),
                    'converged':angle['fits'][method]['converged']})
    axes[0,col].set_title('Axial contrast' if target=='contrast' else 'Central average')
    for ax in axes[:,col]:ax.set_xscale('log');ax.set_xticks([.1,.5,2],['0.1°','0.5°','2°']);ax.set_ylim(-.03,1.05)
    axes[1,col].set_xlabel('Rotation radius (joint pose ball)')
axes[0,0].set_ylabel('Width / no-data width');axes[1,0].set_ylabel('Known-reference correct-sign power')
handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=3,bbox_to_anchor=(.5,1.02),frameon=False,title='EMPIAR geometry; dashed: first order, solid: curvature retained')
fig.tight_layout(rect=[0,0,1,.91]);save(fig,'pose-curvature-width-power')

settings=[('ambient','percentile','#d95f02','Full voxel: bootstrap'),('ambient','same_estimator_ambient_audit','#2166ac','Full voxel: bias audit'),
          ('gaussian123','percentile','#e6ab02','123 modes: bootstrap'),('gaussian123','same_estimator_ambient_audit','#1b9e77','123 modes: bias audit')]
regs=['0.001','0.01','0.1','1.0']
for target in ['center','contrast']:
    fig,axes=plt.subplots(2,3,figsize=(7.2,4.4),sharex='col')
    for col,dataset in enumerate(DATASETS):
        for model,method,color,label in settings:
            data=[]
            for reg in regs:
                source=json.loads((BASE/f'group-bootstrap-reg-{reg}'/f'{dataset}.json').read_text());assert len(source['models'])==3
                block=next(m for m in source['models'] if m['model']==model)
                row=next(r for r in block['results'] if r['target']==target and r['method']==method)
                data.append(row);summary['bootstrap'].append({'dataset':dataset,'regularization':float(reg),'model':model,**row})
            coverage=np.array([r['coverage'] for r in data]);ci=np.array([r['interval_95'] for r in data]).T
            axes[0,col].errorbar(list(map(float,regs)),coverage,yerr=np.abs(ci-coverage),color=color,label=label,marker='o',markersize=3,capsize=2,lw=1)
            axes[1,col].plot(list(map(float,regs)),[r['median_width_fraction_of_no_data'] for r in data],color=color,marker='o',markersize=3,lw=1)
        axes[0,col].set_title('EMPIAR-'+dataset);axes[0,col].axhline(.95,color='black',ls=':',lw=1);axes[0,col].set_ylim(-.03,1.035)
        axes[1,col].set_yscale('log');axes[1,col].set_xlabel('Relative ridge penalty')
        axes[1,col].yaxis.set_major_locator(LogLocator(base=10, numticks=4))
        axes[1,col].yaxis.set_minor_formatter(NullFormatter())
        for ax in axes[:,col]:ax.set_xscale('log');ax.set_xticks([.001,.01,.1,1],['.001','.01','.1','1'])
    axes[0,0].set_ylabel('Coverage');axes[1,0].set_ylabel('Width / no-data width')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=2,bbox_to_anchor=(.5,1.02),frameon=False)
    fig.suptitle(('Central average' if target=='center' else 'Axial contrast')+': fixed poses, independent-map generator',y=1.08)
    fig.tight_layout(rect=[0,0,1,.91]);save(fig,'bootstrap-'+target)

fig,axes=plt.subplots(1,3,figsize=(7.2,2.5))
for ax,dataset in zip(axes,DATASETS):
    source=json.loads((BASE/'reconstruction-comparison'/dataset/'metrics.json').read_text())
    epochs=sorted(map(int,source['epochs']));values=[source['epochs'][str(e)]['prediction']['validation']['nmse'] for e in epochs]
    ax.plot(epochs,values,'o-',color='#2166ac',markersize=3,label='Neural checkpoint')
    for method,color in [('gaussian','#1b9e77'),('voxel','#d95f02')]:
        ax.axhline(source['controls'][method]['validation']['nmse'],ls='--',color=color,label=method.title()+' ridge')
    selected=min(epochs,key=lambda e:(source['epochs'][str(e)]['prediction']['validation']['nmse'],e))
    summary['neural'].append({'dataset':dataset,'available_epochs':epochs,'lowest_tuning_error_epoch':selected,
        'selection_rule':'minimum development tuning NMSE among available checkpoints, ties favor earlier epoch',
        'selected_result':source['epochs'][str(selected)]})
    ax.set_title('EMPIAR-'+dataset);ax.set_xlabel('Training epoch')
axes[0].set_ylabel('Development tuning NMSE');handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=3,frameon=False)
fig.tight_layout(rect=[0,0,1,.9]);save(fig,'neural-convergence')
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(OUT/'README.md').write_text('''# Development comparison figures

Generated by `scripts/summarize_uq_comparisons.py` from complete result files.
All displayed parameter settings are retained. Bootstrap error bars are exact
binomial intervals conditional on the reused resampling-count plan; they are
not voxel-level confidence intervals. The benchmark uses fixed consensus poses,
known Gaussian noise and a deposited-map simulation generator. Neural checkpoint
selection uses development tuning error only. Longer runs may add checkpoints;
this figure is not a proof of global optimization or an end-to-end resolution claim.
''')
