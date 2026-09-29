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
