#!/usr/bin/env python3
"""Publish every prespecified frozen prediction contrast, without reselection."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/confirmation/prediction-v1'
rows=[]; contrasts=[]
for dataset in ['10028','10049','10076']:
    d=json.loads((BASE/dataset/'metrics.json').read_text())
    values=[d['metrics'][m]['nmse'] for m in ['gaussian','voxel','neural']]
    rows.append(dataset+' & '+' & '.join(f'{v:.6f}' for v in values)+r' \\')
    for name,r in d['paired_comparisons'].items():
        lo,hi=r['percentile_bonferroni_six_contrasts']
        contrasts.append(dataset+' & '+name.split('_minus')[0].capitalize()+' & '+f"{r['estimate']:.6f} & [{lo:.6f}, {hi:.6f}]"+r' \\')
text=r'''\begin{table*}[t]
\centering
\caption{Frozen additional-exposure prediction: 4,096 particles per dataset from 115, 71 and 185 source groups absent from development. Models and cohort were locked before image access. Lower normalized prediction MSE is better. Supplied consensus poses/CTFs remain conditional; this does not measure unknown density coverage.}
\label{tab:freshprediction}
\begin{tabular}{lrrr}\toprule
EMPIAR & Gaussian & Voxel & Neural \\\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule\end{tabular}
\medskip

\begin{tabular}{llrl}\toprule
EMPIAR & Comparator minus neural & NMSE difference & Adjusted percentile interval \\\midrule
'''+ '\n'.join(contrasts)+r'''
\bottomrule\end{tabular}
\caption*{All six prespecified paired differences are shown. Intervals use 10,000 whole-exposure bootstrap resamples and a Bonferroni adjustment for six contrasts (99.1667\% marginal percentile intervals). This is a bootstrap approximation, not a finite-sample coverage theorem. No model selection uses these outcomes.}
\end{table*}
'''
# The ICML style does not require caption* support; retain the explanation as
# a small table note instead of adding an incompatible caption package.
text=text.replace(r'\caption*{',r'\par\small{')
(ROOT/'paper/frozen-prediction-results.tex').write_text(text)
