#!/usr/bin/env python3
"""Summarize completed same-weight cross-term audits with source hashes."""
import csv
import hashlib
import json
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def audit_population():
    """List all audit variants/custom folders and the defined source-fit grid."""
    all_attempts=[];unreadable=[]
    for path in sorted(BASE.rglob('*.json')):
        if 'joint-bias-summary' in path.parts:continue
        payload=path.read_bytes()
        try:d=json.loads(payload)
        except json.JSONDecodeError:
            unreadable.append(str(path.relative_to(ROOT)));continue
        if not isinstance(d,dict) or not d.get('stage','').startswith('Exploratory same-weight continuous cross-term audit'):continue
        cfg=d['config'];ok=bool(d.get('complete') and not d.get('error') and not d.get('numerical_failure'))
        all_attempts.append({'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(payload).hexdigest(),
            'source_fit':cfg['fit'],'pose_set':cfg.get('pose_set','joint'),'sharp_cubic':cfg.get('sharp_cubic',False),
            'pilot_pairing':cfg.get('pilot_pairing',False),
            'complete':ok,'error':d.get('error'),'numerical_failure':d.get('numerical_failure',False)})
    population=[]
    for folder in sorted(set(BASE.glob('pose-aware-optimized-shift05*'))|set(BASE.glob('pose-exchange-*'))):
        for path in sorted(folder.glob('*.json')):
            if not re.fullmatch(r'\d{5}-(center|contrast)-[0-9.]+\.json',path.name):continue
            payload=path.read_bytes()
            try:d=json.loads(payload)
            except json.JSONDecodeError:
                unreadable.append(str(path.relative_to(ROOT)));continue
            relative=str(path.relative_to(ROOT));attempts=[a for a in all_attempts if a['source_fit']==relative]
            digest=hashlib.sha256(payload).hexdigest();finished=bool(d.get('complete') and not d.get('error'))
            capture=None
            if not finished:
                # A running fit will change. Preserve the exact inventory input
                # so its current hash does not become an unverifiable pointer.
                snapshot=ROOT/'provenance/uncertainty/result-snapshots'/f'{digest}.json'
                snapshot.parent.mkdir(exist_ok=True)
                if snapshot.exists():assert snapshot.read_bytes()==payload
                else:snapshot.write_bytes(payload)
                capture=str(snapshot.relative_to(ROOT))
            population.append({'source_fit':relative,'sha256':digest,'captured_unfinished_record':capture,
                'fit_complete':finished,'fit_error':d.get('error'),
                'audit_status':'attempted' if attempts else 'unattempted',
                'audits':[a['path'] for a in attempts]})
    return {'all_variants_attempted':len(all_attempts),'all_variants_completed':sum(a['complete'] for a in all_attempts),
        'all_attempts':all_attempts,'source_fit_population':population,
        'population_scope':'All case JSONs in pose-aware-optimized-shift05* and pose-exchange-*; incomplete/unattempted fits explicitly retained',
        'unreadable_json_paths_at_snapshot':sorted(set(unreadable))}


def main():
    out=BASE/'joint-bias-summary';out.mkdir(exist_ok=True);rows=[];sources={};attempts=[]
    for path in sorted((BASE/'joint-bias-audit').glob('*/*.json')):
        d=json.loads(path.read_text())
        sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
        successful=bool(d.get('complete') and not d.get('error') and not d.get('numerical_failure'))
        attempts.append({'path':str(path.relative_to(ROOT)),'complete':successful,
                         'error':d.get('error'),'numerical_failure':d.get('numerical_failure',False)})
        if not successful:continue
        seed=d.get('certificate_seed')
        if seed is None:
            source=ROOT/d['config']['fit']
            if hashlib.sha256(source.read_bytes()).hexdigest()!=d['source_fit_sha256']:
                raise ValueError('Source fit changed; cannot attach historical certificate seed')
            seed=json.loads(source.read_text())['fit']['spectral_upper_certificate']['seed']
        power=min((r['correct_sign_probability'] for r in d['reference_checks']),default=None)
        rows.append({'dataset':d['dataset'],'target':d['target'],'degrees':d['rotation_radius_degrees'],
            'source_fit':path.parent.name,'certificate_seed':seed,'old_width':d['old_selected_relative_half_width'],
            'joint_width':d['new_selected_relative_half_width'],'new_over_old_width':d['relative_width_change'],
            'feasible_over_upper':d.get('existing_feasible_over_new_upper'),'minimum_reference_power':power})
    if not rows:raise RuntimeError('No completed audits')
    refinement=[]
    for path in sorted((BASE/'joint-bias-pilot-sharp-audit').glob('*/*.json')):
        d=json.loads(path.read_text());sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
        if not d.get('complete') or d.get('error') or d.get('numerical_failure'):continue
        sharp_path=BASE/'joint-bias-sharp-audit'/path.parent.name/path.name
        sharp=json.loads(sharp_path.read_text()) if sharp_path.exists() else None
        if sharp is not None:
            assert sharp['complete'] and sharp['source_weights_sha256']==d['source_weights_sha256']
            sources[str(sharp_path.relative_to(ROOT))]=hashlib.sha256(sharp_path.read_bytes()).hexdigest()
        refinement.append({'dataset':d['dataset'],'target':d['target'],'degrees':d['rotation_radius_degrees'],
            'source_fit':path.parent.name,'certificate_seed':d['certificate_seed'],
            'sharp_joint_width':sharp['new_selected_relative_half_width'] if sharp is not None else None,
            'pilot_width':d['new_selected_relative_half_width'],
            'feasible_over_upper':d.get('existing_feasible_over_new_upper'),
            'minimum_reference_power':min((r['correct_sign_probability'] for r in d['reference_checks']),default=None)})
    result={'scope':'Exploratory conditional same-weight post-audits, not density coverage on experimental particles',
            'table_scope':'Default joint-bias-audit folder only; all variants and unattempted source fits listed separately below',
            'completed_cases':len(rows),'attempted_cases':len(attempts),'failed_or_incomplete_cases':len(attempts)-len(rows),
            'attempts':attempts,'certificate_selection':'Every row uses its own source fit certificate; no minimum selected across repeated certificates',
            'pose_set':'Per-particle joint ball: squared scaled rotation norm plus squared scaled translation norm <= 1',
            'source_hashes':sources,'rows':rows,'pilot_refinement_rows':refinement,
            'pilot_refinement_scope':'Separate same-weight, same-event variant using the known pilot and sharp cubic moments; no narrower unknown density class',**audit_population()}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    with (out/'cases.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n');writer.writeheader();writer.writerows(rows)
    lines=[r'\begin{table*}[t]',r'\centering',r'\begin{tabular}{llrrrrr}',r'\toprule',
           r'Stack / fit & Target & Degrees & Old width & Joint width & Stress / upper & Min. power \\',r'\midrule']
    for row in rows:
        fit={'pose-exchange-conic-duals':' / cuts','pose-exchange-adaptive-average':' / full',
             'pose-exchange-sharp-continued':' / sharp fit'}.get(row['source_fit'],'')
        stress='---' if row['feasible_over_upper'] is None else f"{row['feasible_over_upper']:.3f}"
        power=row['minimum_reference_power'];formatted='---' if power is None else (f'{power:.3f}' if power>.001 else f'{power:.1e}')
        if power==0:formatted=r'$0^\dagger$'
        lines.append(f"{row['dataset']}{fit} & {row['target']} & {row['degrees']:g} & {row['old_width']:.3f} & {row['joint_width']:.3f} & {stress} & {formatted} "+r'\\')
    lines.extend([r'\bottomrule',r'\end{tabular}',
        r'\caption{Same-weight residual/pose cross-term audits. Widths are relative to the no-data half-width. All cases use a broad $\sigma_\ell=0.07$ target and a joint scaled pose ball with translation radius $0.5\,\AA$; simultaneous maximum rotation and translation are excluded. Each source cubic remainder and spectral event is retained. Stress ratios reuse saved feasible nonlinear poses; power uses three conditional reference-map scenarios ($\dagger$: numerical underflow). These are exploratory numerical results, not experimental coverage labels.}',
        r'\label{tab:joint-bias}',r'\end{table*}'])
    (ROOT/'paper/tables/joint-bias.tex').write_text('\n'.join(lines)+'\n')
    if refinement:
        with (out/'pilot-cases.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=refinement[0],lineterminator='\n');writer.writeheader();writer.writerows(refinement)
        lines=[r'\begin{table*}[t]',r'\centering',r'\begin{tabular}{llrrrrr}',r'\toprule',
               r'Stack / fit & Target & Degrees & Joint/sharp & Pilot & Stress / upper & Min. power \\',r'\midrule']
        for row in refinement:
            fit={'pose-exchange-conic-duals':' / cuts','pose-exchange-adaptive-average':' / full',
                 'pose-exchange-sharp-continued':' / sharp fit'}.get(row['source_fit'],'')
            sharp='---' if row['sharp_joint_width'] is None else f"{row['sharp_joint_width']:.3f}"
            stress='---' if row['feasible_over_upper'] is None else f"{row['feasible_over_upper']:.3f}"
            power=row['minimum_reference_power'];formatted='---' if power is None else (f'{power:.3f}' if power>.001 else f'{power:.1e}')
            if power==0:formatted=r'$0^\dagger$'
            lines.append(f"{row['dataset']}{fit} & {row['target']} & {row['degrees']:g} & {sharp} & {row['pilot_width']:.3f} & {stress} & {formatted} "+r'\\')
        lines.extend([r'\bottomrule',r'\end{tabular}',
            r'\caption{Known-pilot refinement on the same weights and spectral event as the joint/sharp bound. Both retain the full continuous density class, a $0.07$-field Gaussian width and the joint pose ball with $0.5\,\AA$ translation radius. Widths are relative to the no-data bound; stress and power use the same conditional reference scenarios as Table~\ref{tab:joint-bias}. A width improvement need not yield feature detection ($\dagger$: probability underflow).}',
            r'\label{tab:pilot-pairing}',r'\end{table*}'])
        (ROOT/'paper/tables/pilot-pairing.tex').write_text('\n'.join(lines)+'\n')
    print('DONE',len(rows))


if __name__=='__main__':main()
