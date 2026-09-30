#!/usr/bin/env python3
"""Render every declared CryoLike contrast; no score/setting selection."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    path=ROOT/'results/uncertainty/development/cryolike-comparison-v1/summary.json'
    summary=json.loads(path.read_text())
    if not summary['complete'] or len(summary['records'])!=6:raise ValueError('All six settings required')
    fig,axes=plt.subplots(3,2,figsize=(9.6,8.),layout='constrained')
    methods=['gaussian','voxel','reference'];names=['Gaussian','Voxel','EMDB map']
    rows=[]
    for row,ds in enumerate(['10028','10049','10076']):
        for col,metric in enumerate(['cross_correlation_M','integrated_log_score']):
            ax=axes[row,col]
            for offset,(distance,marker,color) in enumerate([(.4,'o','#0072B2'),(.2,'s','#D55E00')]):
                record=next(r for r in summary['records'] if r['dataset']==ds and r['viewing_distance']==distance)
                if not all(r['scientific_run_complete'] for r in record['cases']):raise ValueError('A method failed')
                for index,method in enumerate(methods):
                    values=record['comparisons'][metric][method+'_minus_neural']
                    estimate=values['estimate'];lo,hi=values['percentile_bonferroni_18']
                    ax.errorbar(estimate,index+(offset-.5)*.2,xerr=[[estimate-lo],[hi-estimate]],
                        marker=marker,color=color,ms=5,capsize=3,lw=1.3,
                        label=f'{distance:.1f} rad' if index==0 else None)
                    rows.append(dict(dataset=ds,metric=metric,viewing_distance=distance,method=method,**values))
            ax.axvline(0,color='#555555',lw=.8,ls='--');ax.set_yticks(range(3),names)
            ax.set_ylim(2.5,-.5);ax.grid(axis='x',alpha=.15)
            ax.set_title(f'EMPIAR-{ds}'+(' (heterogeneous)' if ds=='10076' else ''),fontsize=10)
            ax.set_xlabel('Mean maximum correlation difference' if col==0 else 'Mean integrated log-score difference')
            if row==0:ax.legend(fontsize=8,loc='best')
    fig.suptitle('Original CryoLike scoring: each map minus the locked neural map',fontsize=12)
    folder=ROOT/'paper/figures'
    for suffix in ['pdf','png']:fig.savefig(folder/f'cryolike-comparison.{suffix}',dpi=180)
    plt.close(fig)
    (folder/'cryolike-comparison.json').write_text(json.dumps(dict(complete=True,
        source_sha256=sha(Path(__file__)),summary_sha256=sha(path),contrasts=rows,
        outputs={s:sha(folder/f'cryolike-comparison.{s}') for s in ['pdf','png']},
        interpretation='All 36 paired contrasts, both grids; Bonferroni percentile bootstrap intervals within 18 contrasts per metric. Reused exposures, not density calibration.'),indent=2)+'\n')


if __name__=='__main__':main()
