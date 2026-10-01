#!/usr/bin/env python3
"""Generate the complete six-case folded-width report and paper table."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    path=ROOT/'results/uncertainty/development/folded-ridge-review3-v1/summary.json'
    data=json.loads(path.read_text());assert data['complete'] and len(data['cases'])==6
    rows=data['cases'];assert all('error' not in r for r in rows)
    report=['# Direct folded-width comparison after review 3','',
        'All six declared targets completed in %.2f seconds. All 119 CG solves converged. '
        'The selected widths are .062–.098%% smaller than evaluating the folded critical value '
        'on the previous sum-objective weights. This small improvement does not resolve '
        'pose estimation, density-class calibration, or novelty concerns.'%data['seconds'],'',
        '| Stack | Target | Old width | Selected width | Width ratio | Restricted gap | Global gap | Solves |',
        '|---|---|---:|---:|---:|---:|---:|---:|']
    latex=[r'\begin{table}[t]',r'\centering\small',
        r'\caption{Direct folded-width search at true poses. Ratios compare with the folded width on the previous sum-objective weights. Restricted gaps are below .005; the last column retains the generally looser global diagnostic. No new noisy dataset is drawn.}',
        r'\label{tab:foldedridge}',r'\begin{tabular}{llrrr}',r'\toprule',
        r'Stack & Target & Width & Ratio & Global gap\\',r'\midrule']
    for row in rows:
        f=row['fit'];assert f['all_cg_converged'] and f['converged']
        report.append(f"| {row['dataset']} | {row['target']} | {row['old_sum_weights_folded_width']:.6f} | {f['half_width']:.6f} | {row['width_ratio_to_old']:.6f} | {f['relative_gap']:.6f} | {f['global_relative_gap']:.6f} | {f['evaluations']} |")
        latex.append(f"{row['dataset']} & {row['target']} & {f['half_width']:.4f} & {row['width_ratio_to_old']:.5f} & {f['global_relative_gap']:.4f} \\")
    latex.extend([r'\bottomrule',r'\end{tabular}',r'\end{table}'])
    # Use explicit final backslashes to avoid relying on visually ambiguous escaping.
    for i,line in enumerate(latex):
        if line[:5] in ['10028','10049','10076']:latex[i]=line.rstrip(' \\')+' '+chr(92)*2
    assert sum(row['fit']['evaluations'] for row in rows)==119
    report.extend(['','The search covers [initial ridge / 64, initial ridge * 64]. Its '
        'monotonicity-based relative gap is .00449–.00492 on that finite range. '
        'The global outer-ray diagnostic remains .238 and .289 on 10028, .00492 and '
        '.0332 on 10049, and .00458 and .0571 on 10076. Restricted convergence is '
        'not global optimality. Root evaluation, quadrature and CG errors are '
        'accounted for as described in the protocol; floating-point/NUFFT pads '
        'are diagnostic rather than validated enclosures.','',
        'The known-pose observation design, targets, B=2 and alpha=.05 are '
        'unchanged. No inference images choose the regularization. The source '
        'record retains every ridge evaluation and both matched Gaussian '
        'prior-scale widths. Frozen refitting results are unchanged.',''])
    (ROOT/'research/uncertainty/FOLDED-RIDGE-REVIEW3-RESULTS.md').write_text('\n'.join(report))
    (ROOT/'paper/tables/folded-ridge-review3.tex').write_text('\n'.join(latex)+'\n')


if __name__=='__main__':main()
