#!/usr/bin/env python3
"""Display every paired local-refinement outcome after all frozen attempts finish."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from fourier_splats.uq_end_to_end_summary import METHODS, TEMPLATES, TARGETS
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
DATASETS=['10028','10049','10076']
LABELS=['Fixed folded','Fixed sum','Bounded pose','Mixed 0','Mixed .1','GP 2 fixed','GP 2 pose','GP 1 fixed','GP 1 pose']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    base=ROOT/'results/uncertainty/development/end-to-end-local-pose-summary-v1'
    data={};hashes={}
    for ds in DATASETS:
        p=base/f'{ds}.json';d=json.loads(p.read_text())
        if not d.get('complete') or d['attempted_replicates']!=200 or len(d['groups'])!=72:
            raise ValueError('Every prescribed trial and interval cell required')
        hashes[str(p.relative_to(ROOT))]=sha(p)
        for name,digest in d['input_hashes'].items():
            if sha(ROOT/name)!=digest:
                raise ValueError('Changed frozen record: '+name)
        data[ds]=d
    lookup={(ds,r['template'],r['target'],r['method'],r['image_mode']):r for ds,d in data.items() for r in d['groups']}
    table=[]; report=['# Complete local-pose refitting study','',
        'Exactly 200 independent noisy datasets were attempted for each of three fixed truth/acquisition geometries. All methods, templates, targets and image controls within a dataset are paired. Every failed trial remains in the denominator; interval width summaries describe available results. Exact binomial Monte Carlo intervals and all paired discordance counts are retained in the JSON/CSV summaries.','']
    for ds in DATASETS:
        d=data[ds]
        report += [f"## EMPIAR-{ds}", '', f"Completed {d['completed_replicates']}/200; failures {d['failed_replicates']}. Runtime {d['seconds']:.1f} seconds; peak resident {d['peak_resident_bytes']/1e9:.3f} GB.",'',
            '| Template | Target | Method | Same raw / selected coverage | Independent raw / selected coverage | Median raw / selected relative width | Fallback rate | Independent sign exclusion |',
            '| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |']
        table += [r'\begin{table*}[p]',r'\centering\scriptsize',
            '\\caption{All local-pose outcomes for EMPIAR-'+ds+r'. Coverage is raw/selected, each over 200 planned trials; $W$ is median half-width relative to no data, also raw/selected. Fallback and correct-sign exclusion (independent image, selected) are rates. O/P denote oracle/pilot alignment templates; C/Z denote central/contrast targets. Fixed, bounded and mixed procedures share their fitted weights; GP rows use separate weights. Numerical failures remain in the denominator. Exact binomial Monte Carlo intervals are in the accompanying CSV/JSON.}',
            r'\begin{tabular}{lllrrrrr}',r'\toprule',
            r'Template & Target & Procedure & Same coverage & Independent coverage & $W$ & Fallback & Sign \\',r'\midrule']
        for template in TEMPLATES:
            for target in TARGETS:
                for method,label in zip(METHODS,LABELS):
                    a=lookup[ds,template,target,method,'same_image'];b=lookup[ds,template,target,method,'independent_image']
                    cv=lambda r:f"{r['raw_covered']['fraction']:.3f}/{r['covered']['fraction']:.3f}"
                    wr=a['raw_relative_half_width'].get('median');ws=a['relative_half_width'].get('median')
                    width=f'{wr:.3g}/{ws:.3g}' if wr is not None and ws is not None else 'missing'
                    fb=a['fallback']['fraction'];sign=b['correct_sign_exclusion']['fraction']
                    t='O' if template=='oracle_reference' else 'P';u='C' if target=='center' else 'Z'
                    table.append(f'{t} & {u} & {label} & {cv(a)} & {cv(b)} & {width} & {fb:.3f} & {sign:.3f} '+r'\\')
                    report.append(f'| {t} | {u} | {label} | {cv(a)} | {cv(b)} | {width} | {fb:.3f} | {sign:.3f} |')
                table.append(r'\addlinespace')
        table += [r'\bottomrule',r'\end{tabular}',r'\end{table*}']
    alignment_table=[r'\begin{table*}[t]',r'\centering\footnotesize',
        r'\caption{Complete pose-refinement diagnostics. O/P denote oracle/pilot templates. The score bound is the independent calibration maximum. RMS errors are medians across 200 test datasets. Pose-ball inclusion has an exact binomial 95\% Monte Carlo interval conditional on the realized calibration batch; the rank tolerance theorem instead averages over calibration as well.}',
        r'\begin{tabular}{llrrrr}',r'\toprule',
        r'EMPIAR & Template & Score bound & Rotation RMS ($^\circ$) & Shift RMS (\AA) & Pose-ball inclusion (95\% CI) \\',r'\midrule']
    report += ['', '## Pose calibration and numerical convergence', '',
        '| Stack | Template | Score bound | Median rotation RMS (degrees) | Median shift RMS (A) | Pose-ball inclusion (exact 95% MC interval) |',
        '| --- | --- | ---: | ---: | ---: | --- |']
    for ds,d in data.items():
        for r in d['alignment']:
            t='O' if r['template']=='oracle_reference' else 'P'
            bound=d['calibration_parameters'][r['template']]['joint_score_bound']
            rot=r['rotation_rms_degrees']['median'];shift=r['shift_rms_A']['median']
            p=r['pose_bound_fraction'];lo,hi=r['pose_bound_exact_binomial_95']
            alignment_table.append(f'{ds} & {t} & {bound:.3f} & {rot:.3f} & {shift:.3f} & {p:.3f} ({lo:.3f}, {hi:.3f}) '+r'\\')
            report.append(f'| {ds} | {t} | {bound:.3f} | {rot:.3f} | {shift:.3f} | {p:.3f} ({lo:.3f}, {hi:.3f}) |')
    alignment_table += [r'\bottomrule',r'\end{tabular}',r'\end{table*}']
    table += alignment_table
    for ds,d in data.items():
        report += ['',f"{ds}: {sum(r['converged'] for r in d['numerical'])} converged solves of {sum(r['solves'] for r in d['numerical'])}; every numerical result remains in the coverage evaluation."]
    # Post-outcome arithmetic diagnostic only: decompose every existing bound,
    # without changing a fit, interval, template, or scientific parameter.
    remainder_rows=[]
    report += ['', '## Nonlinear-remainder decomposition', '',
        'Post-outcome reporting diagnostic of every completed fit, with no refitting or selection. Fractions below concern the raw mixed-common-zero width.', '',
        '| Stack | Template | Target | Fits | Median remainder / raw width | Minimum | Maximum |',
        '| --- | --- | --- | ---: | ---: | ---: | ---: |']
    for ds in DATASETS:
        accum={(t,u):[] for t in TEMPLATES for u in TARGETS}
        for rep in range(200):
            p=ROOT/f'results/uncertainty/development/end-to-end-local-pose-v1/{ds}/replicate-{rep:03d}.json'
            rec=json.loads(p.read_text())
            if not rec['complete']:continue
            for t in rec['templates']:
                for f in t['fits']:
                    m=f['mixed_audit']
                    reconstructed=m['density_bias']+m['common_mode_bias']+m['nonlinear_remainder_bias']+m['subgaussian_tail']
                    if not np.isclose(reconstructed,m['half_width'],rtol=1e-12):raise ValueError('Mixed bound decomposition differs')
                    accum[t['template'],f['target']].append(m['nonlinear_remainder_bias']/m['half_width'])
        for (template,target),values in accum.items():
            a=np.asarray(values)
            r=dict(dataset=ds,template=template,target=target,available=len(a),planned=200,
                minimum=float(a.min()) if len(a) else None,median=float(np.median(a)) if len(a) else None,
                maximum=float(a.max()) if len(a) else None)
            remainder_rows.append(r)
            if len(a):report.append(f"| {ds} | {template} | {target} | {len(a)} | {r['median']:.6f} | {r['minimum']:.6f} | {r['maximum']:.6f} |")
    (base/'remainder-diagnostics.json').write_text(json.dumps(dict(complete=True,input_hashes=hashes,
        groups=remainder_rows,source_sha256=sha(Path(__file__)),
        scope='Post-outcome arithmetic decomposition of all existing raw mixed-common-zero bounds; not a new method or predeclared primary outcome.'),indent=2)+'\n')
    (ROOT/'paper/tables/end-to-end-detail.tex').write_text('\n'.join(table)+'\n')
    report += ['', 'A no-data fallback can provide coverage without using inference images or resolving a sign. The deterministic pose construction is simulation-assisted and its marginal tolerance guarantee averages over calibration datasets. Mixed intervals retain unverified conditional centering for estimated designs. These results do not calibrate experimental density coverage or repeat global ab initio reconstruction.']
    (ROOT/'research/uncertainty/END-TO-END-LOCAL-POSE-RESULTS.md').write_text('\n'.join(report)+'\n')
    for template in TEMPLATES:
        # Draw close to the final two-column print width so labels remain
        # legible when the ICML template scales the PDF to textwidth.
        fig,axes=plt.subplots(2,2,figsize=(7.2,6.2),constrained_layout=True)
        configs=[('same_image','raw_covered','Same-image raw coverage',0,1),
                 ('independent_image','raw_covered','Independent-image raw coverage',0,1),
                 ('independent_image','raw_relative_half_width','Median raw width / no data (log scale)',.001,10000),
                 ('independent_image','relative_half_width','Median selected width / no data',0,1.5)]
        for ax,(mode,metric,title,lo,hi) in zip(axes.ravel(),configs):
            matrix=[]
            for ds in DATASETS:
                for target in TARGETS:
                    matrix.append([lookup[ds,template,target,m,mode][metric]['median' if 'half_width' in metric else 'fraction'] for m in METHODS])
            a=np.array(matrix)
            logarithmic=metric=='raw_relative_half_width'
            norm=LogNorm(vmin=lo,vmax=hi) if logarithmic else matplotlib.colors.Normalize(vmin=lo,vmax=hi)
            im=ax.imshow(a,aspect='auto',norm=norm,cmap='viridis')
            ax.set_xticks(range(9),LABELS,rotation=60,ha='right',fontsize=7)
            ax.set_yticks(range(6),[f'{ds} {"C" if t=="center" else "Z"}' for ds in DATASETS for t in TARGETS],fontsize=7)
            ax.set_title(title,fontsize=8)
            for i in range(6):
                for j in range(9):
                    label=f'{a[i,j]:.1e}' if logarithmic else f'{a[i,j]:.2f}'
                    ax.text(j,i,label,ha='center',va='center',fontsize=5.5 if logarithmic else 6.5,
                        color='white' if norm(a[i,j])<.5 else 'black')
            colorbar=fig.colorbar(im,ax=ax,shrink=.7)
            colorbar.ax.tick_params(labelsize=7)
        fig.suptitle(('Oracle-reference' if template=='oracle_reference' else 'Independent-pilot')+' alignment: 200 paired datasets per geometry',fontsize=9)
        for suffix in ['pdf','png']:
            fig.savefig(ROOT/f'paper/figures/end-to-end-{template}.{suffix}',dpi=150)
        plt.close(fig)
    main_table=[r'\begin{figure*}[t]',r'\centering\includegraphics[width=\textwidth]{figures/end-to-end-independent_pilot.pdf}',
        r'\caption{Frozen local-refinement study with independent-pilot alignment. C/Z are central/contrast targets. Every cell uses 200 newly simulated datasets, with poses and weights re-estimated each time. Coverage panels use raw intervals. Width panels distinguish raw bounds (logarithmic color scale) from the no-data-selected procedure. Methods and image controls are paired. Full oracle-template results, selected coverage, binomial Monte Carlo intervals, sign-exclusion and failure counts accompany the figure. These are coarse, fixed-generator simulations, not experimental coverage.}',r'\label{fig:endtoend}',r'\end{figure*}']
    (ROOT/'paper/tables/end-to-end-main.tex').write_text('\n'.join(main_table)+'\n')
    (ROOT/'provenance/uncertainty/end-to-end-reporting-v1.json').write_text(json.dumps(dict(input_hashes=hashes,
        source_snapshot=source_snapshot(ROOT,Path(__file__),[]), complete=True),indent=2)+'\n')
    print('All 600 attempted datasets and 216 procedure cells summarized; raw/fallback outcomes retained',flush=True)


if __name__=='__main__':
    main()
