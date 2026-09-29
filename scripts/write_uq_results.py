#!/usr/bin/env python3
"""Generate draft manuscript tables directly from completed development outputs."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
full=json.loads((ROOT/'results/uncertainty/development/ambient-full-solve.json').read_text())
enriched=json.loads((ROOT/'results/uncertainty/development/ambient-audit.json').read_text())
assert len(full['cases'])==len(enriched['cases'])==3
rows=[];converged=0;seconds=[];enrich_converged=0
for case,other in zip(full['cases'],enriched['cases']):
    assert case['dataset']==other['dataset']
    for width in [.03,.07]:
        targets=[t for t in case['targets'] if t['width_fraction_field']==width]
        e=[t for t in other['targets'] if t['width_fraction_field']==width]
        values=[np.median([t[k] for t in targets]) for k in ['restricted_width_fraction','audited_width_fraction']]
        values.extend([np.median([t['enriched_width_fraction'] for t in e]),np.median([t['matrix_free']['half_width_fraction'] for t in targets])])
        cov=min(r['analytic_coverage'] for t in targets for r in t['coverage'] if r['method']=='restricted_dictionary' and r['truth']=='independent_EMDB_voxel_reference')
        rows.append(case['dataset']+' & '+f'{width:.2f}'+' & '+' & '.join(f'{v:.3f}' for v in values)+r' \\')
        for t in targets:
            converged+=t['matrix_free']['converged'];seconds.append(t['matrix_free']['seconds'])
        enrich_converged+=sum(t['converged'] for t in e)
        print(case['dataset'],width,values,'minimum dictionary/reference coverage',cov)
text=r'''\begin{table*}[t]
\centering
\caption{Fixed-pose ambient audit: median interval width relative to the no-data width across two density targets. All three methods after auditing use the same supported voxel class. The first column instead assumes the initial 33-dimensional dictionary is correct and is not an honest interval for that larger class. All entries are development results.}
\label{tab:ambient}
\begin{tabular}{llrrrr}
\toprule
Dataset & Averaging width / field & Dictionary only & Same weights, audited & Enriched & Full matrix-free \\
\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule
\end{tabular}
\end{table*}

'''
text+=f'Table~\\ref{{tab:ambient}} reports median widths across the two targets. The full-space solver reached its 0.5\\% objective-gap tolerance in {converged} of 12 cases, taking {min(seconds):.1f}--{max(seconds):.1f} seconds per target on the Mac. The enrichment solver reached the same outer tolerance in {enrich_converged} of 12 cases within 100 rounds; nonconverged results retain valid bounds but do not establish optimality. The solvers optimize the conservative sum-width, so their sharper reported widths need not be ordered exactly by their objective gaps.\n\n'
text+=r'''Dictionary-only intervals can be narrow because they exclude admissible density errors. For the EMPIAR-10028 central broad average, analytic coverage on the independent voxel-map generator is approximately $4.2\times10^{-5}$. The same image weights, with the ambient bias bound, cover at the prescribed level over the ambient ball. On each estimator's constructed ambient bias boundary, the audited intervals attain nominal 0.95 Gaussian coverage up to numerical tolerance. These boundary calculations check the theorem and implementation; they are not independent evidence of practical class calibration. Broader targets have substantially smaller relative widths than fine targets. The prescribed norm radius and fixed-pose assumption are essential to interpreting that result.
'''
(ROOT/'paper/uncertainty-results.tex').write_text(text)

pose=json.loads((ROOT/'results/uncertainty/development/ambient-pose.json').read_text())
refined=json.loads((ROOT/'results/uncertainty/development/ambient-pose-large-radius-refinement.json').read_text())
refinement={c['dataset']:c['angles'][0]['fits']['bounded_pose_ambient'] for c in refined['cases']}
pose_rows=[]
for case in pose['cases']:
    values=[a['fits']['bounded_pose_ambient']['width_fraction_of_no_data'] for a in case['angles']]
    final=refinement[case['dataset']]
    pose_rows.append(case['dataset']+' & '+' & '.join(f'{v:.3f}' for v in values)+f" & {final['width_fraction_of_no_data']:.3f}"+r' \\')
pose_text=r'''\begin{table}[t]
\centering
\caption{Joint ambient/pose interval width divided by the no-data width. Columns 2--4 use a 100-iteration cap. The last column repeats $2^\circ$ with up to 500 iterations; EMPIAR-10049 retains a 0.548\% gap, above the 0.5\% tolerance. Others meet tolerance. These are different optimization budgets, not separate independent replications.}
\label{tab:ambientpose}
\begin{tabular}{lrrrr}
\toprule
EMPIAR & $0.1^\circ$ & $0.5^\circ$ & $2^\circ$ & $2^\circ$, longer \\
\midrule
'''+ '\n'.join(pose_rows)+r'''
\bottomrule
\end{tabular}
\end{table}
'''
(ROOT/'paper/uncertainty-pose-results.tex').write_text(pose_text)

power_rows=[]
for dataset in ['10028','10049','10076']:
    values=[]
    for filename in ['ambient-pose-power.json','ambient-pose-center-power.json']:
        source=json.loads((ROOT/'results/uncertainty/development'/filename).read_text())
        case=next(c for c in source['cases'] if c['dataset']==dataset)
        for angle in case['angles']:
            r=next(r for r in angle['coverage'] if r['method']=='bounded_pose_ambient' and r['truth']=='independent_EMDB_reference_coherent_pose')
            values.append(r['correct_sign_probability'])
    power_rows.append(dataset+' & '+' & '.join(f'{v:.3f}' for v in values)+r' \\')
(ROOT/'paper/uncertainty-power-results.tex').write_text(r'''\begin{table}[t]
\centering
\caption{Probability that the nonlinear interval establishes the correct sign under the independent deposited-map generator and a coherent allowed rotation. Targets are broad averages or contrasts; these are not atomic-resolution claims. High coverage alone does not ensure detection.}
\label{tab:power}
\begin{tabular}{lrrrr}
\toprule
 & \multicolumn{2}{c}{Axial contrast} & \multicolumn{2}{c}{Central average} \\
EMPIAR & $0.1^\circ$ & $0.5^\circ$ & $0.1^\circ$ & $0.5^\circ$ \\
\midrule
'''+ '\n'.join(power_rows)+r'''
\bottomrule
\end{tabular}
\end{table}
''')

reconstruction_rows=[]
for dataset in ['10028','10049','10076']:
    source=json.loads((ROOT/'results/uncertainty/development/reconstruction-comparison'/dataset/'metrics.json').read_text())
    values=[source['controls'][method]['test']['nmse'] for method in ['gaussian','voxel']]
    values.append(source['epochs']['20']['prediction']['test']['nmse'])
    limit=source['epochs']['20']['conditional_half_fsc_resolution']['angstrom']
    reconstruction_rows.append(dataset+' & '+' & '.join(f'{v:.4f}' for v in values)+f' & {limit:.2f}'+r' \\')
(ROOT/'paper/uncertainty-reconstruction-results.tex').write_text(r'''\begin{table}[t]
\centering
\caption{Experimental-image prediction error on development test groups (normalized MSE; lower is better). Neural runs are 20-epoch checkpoints; convergence is still being assessed. All half-map FSCs remain above 0.143 through the sampled limit in the last column, which is not a measured crossing.}
\label{tab:reconstruction}
\begin{tabular}{lrrrr}
\toprule
EMPIAR & Gaussian & Voxel & Neural & Limit (\AA) \\
\midrule
'''+ '\n'.join(reconstruction_rows)+r'''
\bottomrule
\end{tabular}
\end{table}
''')
