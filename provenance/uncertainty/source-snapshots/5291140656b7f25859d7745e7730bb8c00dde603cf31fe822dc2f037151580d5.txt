#!/usr/bin/env python3
"""All-case matched-prior report; incomplete fits and hash changes are explicit."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
DATASETS = ['10028', '10049', '10076']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def span(x):
    a = np.asarray(x, float)
    if not len(a) or not np.isfinite(a).all():
        raise ValueError('Nonempty finite summary required')
    return dict(minimum=float(a.min()), median=float(np.median(a)), maximum=float(a.max()))


def main():
    directory = BASE/'continuous-gaussian-review2-v2'
    records, hashes, source = [], {}, {}
    for ds in DATASETS:
        p = directory/f'{ds}.json'; d = json.loads(p.read_text())
        if not d.get('complete') or len(d['records']) != 16:
            raise ValueError(f'{ds}: require all 16 prescribed fits')
        hashes[str(p.relative_to(ROOT))] = sha(p); source[ds] = d
        for name, expected in d['input_hashes'].items():
            if sha(ROOT/name) != expected:
                raise ValueError(f'Changed input: {name}')
            hashes[name] = expected
        ap = directory/f'{ds}-design.npz'
        if sha(ap) != d['design_arrays_sha256']:
            raise ValueError('Changed design')
        hashes[str(ap.relative_to(ROOT))] = sha(ap)
        expected = {(t, tau, pose) for t in ['pilot_region_1', 'pilot_region_2', 'pilot_region_3', 'matched_center']
                    for tau in [1., 2.] for pose in ['fixed_pose', 'local_gaussian_pose']}
        seen = set()
        for row in d['records']:
            key = row['target'], row['prior_directional_sd'], row['pose_model']
            if key not in expected or key in seen or len(row['checks']) != 14:
                raise ValueError('Missing or duplicated fit/scenario')
            seen.add(key)
            cases = {(frame, setting) for frame in ['original', 'registered'] for setting in
                ['nominal']+[f'{angle}deg_{direction}' for angle in [0, 1, 2] for direction in ['coherent_x', 'random_boundary']]}
            if {(r['frame'], r['scenario']) for r in row['checks']} != cases:
                raise ValueError('Missing declared frame/scenario')
            wp = directory/f'{ds}-{key[0]}-tau{key[1]:g}-{key[2]}.npz'
            if sha(wp) != row['weights_sha256']:
                raise ValueError('Changed weights')
            hashes[str(wp.relative_to(ROOT))] = sha(wp)
            for check in row['checks']:
                records.append(dict(dataset=ds, target=key[0], tau=key[1], pose_model=key[2],
                    converged=row['fit']['converged'],
                    cg_relative_residual=row['fit']['cg_relative_residual'],
                    cg_iterations=row['fit']['cg_iterations'],
                    variance_identity_relative_error=row['fit']['variance_identity_relative_error'],
                    width_over_fixed_audit=row['half_width_over_fixed_audit'], **check))
        if seen != expected:
            raise ValueError('Missing declared fit')
    output = BASE/'continuous-gaussian-review2-summary-v2'
    if output.exists():
        raise RuntimeError('Preserve prior final report')
    output.mkdir(parents=True)
    summary = dict(complete=True, fits=48, conditional_scenarios=len(records), input_hashes=hashes,
        source_snapshot=source_snapshot(ROOT, Path(__file__), []),
        scope='Conditional prescribed-noise calculations, not refitted-pose or experimental coverage. No failed fit omitted.', groups=[],
        resources={ds:{k:source[ds][k] for k in ['seconds','setup_seconds','peak_resident_bytes','preconditioner_diagnostics']} for ds in DATASETS})
    for ds in DATASETS:
        for tau in [2., 1.]:
            for pose in ['fixed_pose', 'local_gaussian_pose']:
                rr = [r for r in records if (r['dataset'],r['tau'],r['pose_model']) == (ds,tau,pose)]
                nominal = [r for r in rr if r['frame']=='registered' and r['scenario']=='nominal']
                group = dict(dataset=ds, tau=tau, pose_model=pose, fits=4,
                    converged_fits=sum(r['converged'] for r in nominal),
                    width_over_fixed_audit=span([r['width_over_fixed_audit'] for r in nominal]),
                    half_width_over_no_data=span([r['half_width_over_no_data'] for r in nominal]),
                    cg_relative_residual=span([r['cg_relative_residual'] for r in nominal]),
                    cg_iterations=span([r['cg_iterations'] for r in nominal]),
                    registered_nominal_coverage=span([r['analytic_fixed_signal_coverage'] for r in nominal]),
                    registered_nominal_sign_probability=span([r['correct_sign_probability'] for r in nominal]),
                    both_frames_all_scenarios_coverage=span([r['analytic_fixed_signal_coverage'] for r in rr]),
                    both_frames_all_scenarios_sign_probability=span([r['correct_sign_probability'] for r in rr]))
                summary['groups'].append(group)
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (output/'all-cases.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    table = [r'\begin{table}[t]',r'\centering\small',
        r'\caption{Matched continuous Gaussian comparison: range over four locked features per stack. $R$ is half-width divided by the original fixed-pose audit width. $C_{\min}$ is minimum analytic coverage across both generator frames and all seven prescribed pose settings. All 48 solves converge; this is conditional supplied-noise coverage, not refitted-pose validation.}',
        r'\label{tab:matchedgaussian}',r'\begin{tabular}{llrrr}',r'\toprule',
        r'Stack & Pose model & $\tau$ & $R$ range & $C_{\min}$ \\',r'\midrule']
    report = ['# Complete matched continuous Gaussian comparison','',
        'All 48 prescribed fits and 672 conditional frame/scenario cases are retained. Source, design and weight hashes were verified. These calculations hold poses/weights fixed; they are not experimental density-coverage evidence.','',
        '| Stack | Pose covariance | Prior scale | Width / audit range | Min coverage, all cases | Max nominal registered sign probability |',
        '| --- | --- | ---: | ---: | ---: | ---: |']
    for g in summary['groups']:
        r=g['width_over_fixed_audit'];c=g['both_frames_all_scenarios_coverage']['minimum'];p=g['registered_nominal_sign_probability']['maximum']
        pose='Fixed' if g['pose_model']=='fixed_pose' else 'Linearized'
        table.append(f"{g['dataset']} & {pose} & {g['tau']:g} & {r['minimum']:.3f}--{r['maximum']:.3f} & {c:.3f} "+r'\\')
        report.append(f"| {g['dataset']} | {pose} | {g['tau']:g} | {r['minimum']:.6g}–{r['maximum']:.6g} | {c:.8g} | {p:.6g} |")
    table += [r'\bottomrule',r'\end{tabular}',r'\end{table}']
    # Do not silently print a convergence claim if any future input differs.
    if not all(r['converged'] for r in records):
        table[2]=table[2].replace('All 48 solves converge;', 'Some solves did not converge;')
    (ROOT/'paper/tables/continuous-gaussian-v2-main.tex').write_text('\n'.join(table)+'\n')
    report += ['', 'The scale-matched prior can share the fixed-pose uniform guarantee by the proved normal-envelope inequality; the half-scale prior need not. Favorable coverage on the two specific generators does not validate the class or a nonlinear Gaussian-pose approximation.',
        '', 'The interrupted rank-1,024 run remains incomplete, with all four nonconverged fits and the numerical amendment preserved. The complete rank-8,192 rerun changes preconditioning and transform caching only.']
    (ROOT/'research/uncertainty/CONTINUOUS-GAUSSIAN-V2-RESULTS.md').write_text('\n'.join(report)+'\n')
    fig, axes=plt.subplots(1,3,figsize=(12,3.7),sharey=True,constrained_layout=True)
    colors=['#3569a8','#91b5d3','#d66b36','#ebaf7e']
    methods=[(2.,'fixed_pose'),(2.,'local_gaussian_pose'),(1.,'fixed_pose'),(1.,'local_gaussian_pose')]
    for ax, ds in zip(axes,DATASETS):
        for j,(tau,pose) in enumerate(methods):
            rr=[r for r in records if (r['dataset'],r['tau'],r['pose_model'],r['frame'],r['scenario'])==(ds,tau,pose,'registered','nominal')]
            ax.scatter(np.arange(4)+(j-1.5)*.13,[r['width_over_fixed_audit'] for r in rr],s=27,color=colors[j],label=f'τ={tau:g}, '+('fixed' if pose=='fixed_pose' else 'linearized pose'))
        ax.axhline(1,color='.4',lw=1,ls='--');ax.set_title('EMPIAR-'+ds);ax.set_xticks(range(4),['Region 1','Region 2','Region 3','Center'],rotation=20)
    axes[0].set_ylabel('Gaussian half-width / fixed-pose audit half-width');axes[-1].legend(fontsize=7)
    for suffix in ['pdf','png']:
        fig.savefig(ROOT/f'paper/figures/continuous-gaussian-v2.{suffix}',dpi=150)
    plt.close(fig)
    print('Verified',len(hashes),'inputs;',len(records),'conditional cases;',sum(g['converged_fits'] for g in summary['groups']),'converged fits',flush=True)


if __name__ == '__main__':
    main()
