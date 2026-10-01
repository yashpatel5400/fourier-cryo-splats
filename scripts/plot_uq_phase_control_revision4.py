#!/usr/bin/env python3
"""Presentation-only redraw of the unchanged phase-control results."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'results/uncertainty/development/phase-split-control-v1/summary.json').read_text())
labels=['Known phase','Same-image aligned','Independent aligned: nominal SD','Independent aligned: sample SD','Invariant cross-power bound']
colors=['#555555','#d55e00','#0072b2','#56b4e9','#009e73']
plt.rcParams.update({'pdf.fonttype':42})
fig,axes=plt.subplots(2,4,figsize=(9,4.9),sharex=True)
for col,sigma in enumerate([2.,1.,.5,.25]):
 cases=[r for r in d['records'] if r['sigma']==sigma];ns=[r['particles'] for r in cases]
 for mi,(label,color) in enumerate(zip(labels,colors)):
  records=[r['methods'][mi] for r in cases];means=np.array([r['coverage'] for r in records]);bounds=np.array([r['exact_binomial_95'] for r in records])
  axes[0,col].errorbar(ns,means,yerr=np.vstack([means-bounds[:,0],bounds[:,1]-means]),label=label,color=color,marker='o',ms=3,lw=1,capsize=1)
  axes[1,col].plot(ns,[r['median_half_width'] for r in records],color=color,marker='o',ms=3,lw=1)
 axes[0,col].axhline(.95,color='black',linestyle=':',lw=.8);axes[0,col].set_title(f'Amplitude / noise = {1/sigma:g}',fontsize=9);axes[0,col].set_ylim(-.035,1.04)
 axes[1,col].set_yscale('log');axes[1,col].yaxis.set_major_locator(LogLocator(base=10,numticks=4));axes[1,col].yaxis.set_minor_formatter(NullFormatter())
 axes[1,col].set_xlabel('Independent particles',fontsize=8)
 for ax in axes[:,col]:
  ax.set_xscale('log',base=4);ax.set_xticks([16,256,4096],['16','256','4096']);ax.tick_params(labelsize=7);ax.grid(alpha=.2)
axes[0,0].set_ylabel('Coverage of true amplitude',fontsize=8);axes[1,0].set_ylabel('Median half-width',fontsize=8)
h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,loc='lower center',ncol=2,fontsize=8,frameon=False)
fig.tight_layout(rect=[0,.17,1,1]);fig.savefig(ROOT/'paper/figures/phase-split-control.pdf');fig.savefig(ROOT/'paper/figures/phase-split-control.png',dpi=160)
