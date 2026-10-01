#!/usr/bin/env python3
"""Summarize every declared acquisition diagnostic; retain the failed transfer."""
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from review_uq_candidate import snapshot_references

ROOT=Path(__file__).resolve().parents[1]


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    folder=ROOT/'results/uncertainty/development/raw-movie-pilot-v1'
    p=folder/'summary.json';d=json.loads(p.read_text())
    if not d.get('complete') or len(d['frames'])!=16 or len(d['patch_differences'])!=512 or len(d['disjoint_pair_correlations'])!=1792:
        raise ValueError('Complete predeclared acquisition analysis required')
    hashes={str(p.relative_to(ROOT)):sha(p),**d['input_hashes']}
    hashes[str((folder/'diagnostics.npz').relative_to(ROOT))]=d['arrays_sha256']
    hashes.update({str((folder/name).relative_to(ROOT)):h for name,h in d['figure_hashes'].items()})
    for name,h in hashes.items():
        if sha(ROOT/name)!=h:raise ValueError('Changed acquisition input/output: '+name)
    canonical=ROOT/'provenance/uncertainty/raw-movie-pilot-v1.json'
    source=json.loads(canonical.read_text())
    for key in ['previous_failed_attempt','resume_record']:
        if sha(ROOT/source[key])!=source[key+'_sha256']:raise ValueError('Changed recovery evidence')
        hashes[source[key]]=source[key+'_sha256']
    levels=[0,.05,.5,.95,1]
    rows={key:np.quantile([r[key] for r in d['patch_differences']],levels).tolist()
        for key in ['sd','horizontal_correlation','vertical_correlation']}
    rows['disjoint_difference_correlation']=np.quantile([r['correlation'] for r in d['disjoint_pair_correlations']],levels).tolist()
    out=folder/'report.json'
    if out.exists():raise RuntimeError('Preserve existing report')
    record=dict(complete=True,input_hashes=hashes,quantile_levels=levels,quantiles=rows,
        source_sha256=sha(Path(__file__)),scope=d['scope'],
        frame_mean_range=[min(r['mean'] for r in d['frames']),max(r['mean'] for r in d['frames'])],
        frame_sd_range=[min(r['sd'] for r in d['frames']),max(r['sd'] for r in d['frames'])],
        maximum_parseval_error=d['maximum_parseval_error'],raw_movie_sha256=source['sha256'])
    out.write_text(json.dumps(record,indent=2)+'\n')
    md=['# Complete one-movie acquisition pilot','',
        'One preselected EMPIAR-10028 movie contains 16 frames of 4096 x 4096 float32 pixels. The first transfer broke after 922,746,880 saved bytes. A range request recovered the remaining 150,995,968 bytes after checking ETag, Last-Modified, Content-Range and size. The failed record, original partial file and recovery tail remain preserved. The assembled 1,073,742,848-byte movie has SHA-256 `'+source['sha256']+'`.','',
        'All 64 fixed patches and eight disjoint frame pairs were analyzed: 512 patch differences and 1,792 within-patch cross-pair correlations. No patches or frames were selected by outcome.','',
        '| Statistic | Minimum | 5th percentile | Median | 95th percentile | Maximum |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for key,values in rows.items():md.append('| '+key.replace('_',' ')+' | '+' | '.join(f'{x:.6g}' for x in values)+' |')
    md += ['',f"Frame means range from {record['frame_mean_range'][0]:.3f} to {record['frame_mean_range'][1]:.3f}; frame SDs range from {record['frame_sd_range'][0]:.3f} to {record['frame_sd_range'][1]:.3f}, in deposited pixel units. Maximum absolute Parseval discrepancy is {d['maximum_parseval_error']:.3g}.",'',
        'Spatial correlation is substantial. Small correlations between disjoint differences are descriptive and do not prove independence, Gaussianity, or covariance transfer. These differences retain dose, motion and detector contributions. The single movie does not validate a homogeneous density, calibrate poses, or resolve the source-group identity limitation. No density intervals or reconstruction are derived from this pilot.','',
        'The original movie remains at its recorded EMPIAR URL rather than being duplicated in the release. The release includes every diagnostic, display, protocol, transfer record and analysis source snapshot. The original failed prefix and tail are transport artifacts, not extra scientific observations.']
    (ROOT/'research/uncertainty/RAW-MOVIE-PILOT-RESULTS.md').write_text('\n'.join(md)+'\n')
    tex=[r'\begin{table*}[t]',r'\centering\small',
        r'\caption{Predeclared acquisition diagnostics from one raw 10028 movie. Difference statistics use every one of 512 patch/pair combinations; disjoint-pair correlations use all 1,792 comparisons. Quantiles are descriptive across dependent summaries, not independent-replicate confidence intervals.}',
        r'\begin{tabular}{lrrrrr}',r'\toprule',r'Statistic & Minimum & 5th percentile & Median & 95th percentile & Maximum \\',r'\midrule']
    labels=['Difference SD','Horizontal neighbor correlation','Vertical neighbor correlation','Disjoint-difference correlation']
    for label,values in zip(labels,rows.values()):tex.append(label+' & '+' & '.join(f'{x:.4g}' for x in values)+r' \\')
    tex += [r'\bottomrule',r'\end{tabular}',r'\end{table*}']
    (ROOT/'paper/tables/raw-movie-diagnostics.tex').write_text('\n'.join(tex)+'\n')
    # Package analysis products and exact provenance, not a second copy of the
    # publicly hosted 1 GiB movie or incomplete transport bytes.
    members={p,out,*folder.glob('*.npz'),*folder.glob('*.pdf'),*folder.glob('*.png'),canonical,
        ROOT/source['previous_failed_attempt'],ROOT/source['resume_record'],
        ROOT/'research/uncertainty/RAW-MOVIE-PILOT-PROTOCOL.md',
        ROOT/'research/uncertainty/RAW-MOVIE-PILOT-RESULTS.md',Path(__file__)}
    for q in list(members):
        if q.suffix=='.json':
            for ref in snapshot_references(json.loads(q.read_text())):
                s=ROOT/ref
                if sha(s)!=s.stem:raise ValueError('Changed archived analysis source')
                members.add(s)
    owner=str(p.relative_to(ROOT));oh=sha(p)
    spec=dict(version='v0.7.0-dev-movie-diagnostics',files=[dict(path=str(q.relative_to(ROOT)),
        sha256=sha(q),owner=owner,owner_sha256=oh) for q in sorted(members)],
        required_prior_arrays=[{'url':source['url'],'bytes':source['expected_bytes'],'sha256':source['sha256'],'path':source['path']}],
        scope=d['scope'],reproduction='Download the single source movie and verify its recorded size/header/hash. Use the archived diagnostic source. Preserve an existing output directory before an independent rerun; the runner intentionally refuses overwrites. No particle picking or pose processing is part of this analysis.')
    (ROOT/'provenance/uncertainty/raw-movie-artifact-spec-v1.json').write_text(json.dumps(spec,indent=2)+'\n')
    print('Verified 16 frames, 512 differences and 1792 correlations;',len(members),'artifact members')


if __name__=='__main__':main()
