#!/usr/bin/env python3
"""Report the full completed target family and its frozen noise recalibration."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--report-version', type=int, default=1)
    args = parser.parse_args()
    if args.report_version < 1: raise ValueError('Positive report version required')
    old_path = ROOT/'results/uncertainty/development/pilot-selected-summary-v1/summary.json'
    fresh_dir = ROOT/'results/uncertainty/confirmation/noise-calibration-v1'
    fresh_path = fresh_dir/'summary.json'
    old = json.loads(old_path.read_text()); fresh = json.loads(fresh_path.read_text())
    assert old['complete'] and fresh['complete']
    assert len(old['conditional']) == len(old['experimental']) == len(fresh['records']) == 48
    hashes = {str(p.relative_to(ROOT)):sha(p) for p in [old_path,fresh_path]}
    for name, digest in old['source_hashes'].items():
        assert sha(ROOT/name) == digest, name
    cases = {}
    for ds in ['10028','10049','10076']:
        p = fresh_dir/f'{ds}.json'; assert sha(p) == fresh['source_hashes'][p.name]
        cases[ds] = json.loads(p.read_text()); assert cases[ds]['complete'] and not cases[ds].get('error')
        hashes[str(p.relative_to(ROOT))] = sha(p)
    out = ROOT/f'results/uncertainty/development/locked-target-report-v{args.report_version}'
    if out.exists(): raise RuntimeError('Preserve earlier report, including failures')
    out.mkdir(parents=True)
    poses = ['fixed_pose','shift_only','pose_1deg','pose_2deg']
    rows = []
    for ds in cases:
        for pose in poses:
            c = [r for r in old['conditional'] if r['dataset']==ds and r['pose_class']==pose]
            e = [r for r in fresh['records'] if r['dataset']==ds and r['pose_class']==pose]
            assert len(c)==len(e)==4
            rows.append({'dataset':ds,'pose_class':pose,
                'known_noise_width_range':[min(r['relative_width'] for r in c),max(r['relative_width'] for r in c)],
                'known_noise_reference_power_range':[min(r['minimum_reference_power'] for r in c),max(r['minimum_reference_power'] for r in c)],
                'fresh_calibrated_width_range':[min(r['relative_half_width'] for r in e),max(r['relative_half_width'] for r in e)],
                'zero_exclusion_features':[r['feature'] for r in e if r['excludes_zero']],
                'fallback_count':sum(r['uses_no_data'] for r in e),
                'approximate_reference_inside_count':sum(r['approximate_reference_inside'] for r in e)})
    result = {'complete':True,'source_hashes':hashes,'rows':rows,
        'fresh_to_old_sd_ranges':{ds:[min(f['fresh_to_old_sd_ratio'] for f in d['features']),
                                    max(f['fresh_to_old_sd_ratio'] for f in d['features'])] for ds,d in cases.items()},
        'scope':'Complete descriptive summary; no new estimation, selection or calibration. Each pose class has its own 12-feature family budget.',
        'reference_note':'Approximate-map inclusion is not empirical density coverage; class/noise/pose assumptions remain unverified.'}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    # Show every feature and pose class, with the same normalization as the method.
    plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral'],'pdf.fonttype':42})
    fig, axes = plt.subplots(1,3,figsize=(7.2,2.85),sharey=True)
    colors = ['#0072B2','#009E73','#E69F00','#CC79A7']
    names = ['Fixed poses','0.5 Å shifts','1° + 0.5 Å','2° + 0.5 Å']
    for ax,(ds,case) in zip(axes,cases.items()):
        for i,feature in enumerate(case['features']):
            for j,r in enumerate(feature['records']):
                assert r['pose_class']==poses[j]
                original = next(x for x in old['experimental'] if x['dataset']==ds and x['feature']==feature['target'] and x['pose_class']==poses[j])
                denominator = original['no_data_half_width']
                ax.errorbar(r['interval_center']/denominator,i+(j-1.5)*.15,
                    xerr=r['interval_half_width']/denominator,fmt='o',color=colors[j],
                    markersize=2.5,capsize=1.6,linewidth=.85,label=names[j] if i==0 else None)
            ref = feature['records'][0]['approximate_reference_target']/denominator
            ax.plot([ref,ref],[i-.32,i+.32],'k:',linewidth=.9,label='Approx. map' if i==0 else None)
        ax.axvline(0,color='.55',linewidth=.5);ax.grid(axis='x',alpha=.18)
        ax.set_title(f'EMPIAR-{ds}'+(' (heterogeneous)' if ds=='10076' else ''),fontsize=9)
        ax.set_xlabel('Feature / no-data half-width',fontsize=8)
        ax.tick_params(labelsize=8)
    axes[0].set_yticks(range(4),['Region 1','Region 2','Region 3','Center']);axes[0].invert_yaxis()
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=5,fontsize=7,frameon=False)
    fig.tight_layout(rect=(0,.10,1,1))
    for ext in ['pdf','png']:
        fig.savefig(ROOT/f'paper/figures/fresh-noise-targets.{ext}',dpi=200)
    plt.close(fig)
    table = [r'\begin{table*}[t]',r'\centering\small',
        r'\caption{Frozen fresh-exposure noise recalibration of twelve fixed estimators. Each entry gives the range of relative half-widths across four features and the number excluding zero. These are conditional sensitivity results, not calibrated experimental density coverage. All approximate map values lie inside the intervals.}',
        r'\label{tab:fresh-targets}',r'\begin{tabular}{lcccc}',r'\toprule',
        r'Stack & Fixed poses & 0.5 \AA\ shifts & $1^\circ+0.5$ \AA & $2^\circ+0.5$ \AA \\',r'\midrule']
    for ds in cases:
        values=[]
        for r in [r for r in rows if r['dataset']==ds]:
            a,b=r['fresh_calibrated_width_range'];interval=f'{a:.3f}' if f'{a:.3f}'==f'{b:.3f}' else f'{a:.3f}--{b:.3f}'
            values.append(f'{interval} ({len(r["zero_exclusion_features"])}/4)')
        table.append(ds+' & '+' & '.join(values)+r' \\')
    table += [r'\bottomrule',r'\end{tabular}',r'\end{table*}']
    (ROOT/'paper/tables/fresh-noise-targets.tex').write_text('\n'.join(table)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__': main()
