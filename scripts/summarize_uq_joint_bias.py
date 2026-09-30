#!/usr/bin/env python3
"""Summarize completed same-weight cross-term audits with source hashes."""
import csv
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    out=BASE/'joint-bias-summary';out.mkdir(exist_ok=True);rows=[];sources={}
    for path in sorted((BASE/'joint-bias-audit').glob('*/*.json')):
        d=json.loads(path.read_text())
        if not d.get('complete'):continue
        sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
        power=min((r['correct_sign_probability'] for r in d['reference_checks']),default=None)
        rows.append({'dataset':d['dataset'],'target':d['target'],'degrees':d['rotation_radius_degrees'],
            'source_fit':path.parent.name,'old_width':d['old_selected_relative_half_width'],
            'joint_width':d['new_selected_relative_half_width'],'new_over_old_width':d['relative_width_change'],
            'feasible_over_upper':d.get('existing_feasible_over_new_upper'),'minimum_reference_power':power})
    if not rows:raise RuntimeError('No completed audits')
    result={'scope':'Exploratory conditional same-weight post-audits, not density coverage on experimental particles',
            'completed_cases':len(rows),'source_hashes':sources,'rows':rows}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    with (out/'cases.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n');writer.writeheader();writer.writerows(rows)
    lines=[r'\begin{table*}[t]',r'\centering',r'\begin{tabular}{llrrrrr}',r'\toprule',
           r'Stack / fit & Target & Degrees & Old width & Joint width & Stress / upper & Min. power \\',r'\midrule']
    for row in rows:
        fit=' / cuts' if row['source_fit']=='pose-exchange-conic-duals' else ''
        stress='---' if row['feasible_over_upper'] is None else f"{row['feasible_over_upper']:.3f}"
        power=row['minimum_reference_power'];formatted='---' if power is None else (f'{power:.3f}' if power>.001 else f'{power:.1e}')
        lines.append(f"{row['dataset']}{fit} & {row['target']} & {row['degrees']:g} & {row['old_width']:.3f} & {row['joint_width']:.3f} & {stress} & {formatted} "+r'\\')
    lines.extend([r'\bottomrule',r'\end{tabular}',
        r'\caption{Same-weight residual/pose cross-term audits. Widths are relative to the no-data half-width. All cases use a broad $\sigma_\ell=0.07$ target and $0.5\,\AA$ translations. The old cubic remainder and original spectral event are retained here. Stress ratios reuse saved feasible nonlinear poses; power uses the three conditional reference-map scenarios. These are exploratory numerical results, not experimental coverage labels.}',
        r'\label{tab:joint-bias}',r'\end{table*}'])
    (ROOT/'paper/tables/joint-bias.tex').write_text('\n'.join(lines)+'\n')
    print('DONE',len(rows))


if __name__=='__main__':main()
