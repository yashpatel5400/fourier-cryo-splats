#!/usr/bin/env python3
"""Package the complete conditional-noise/preferred-view increment."""
import hashlib,json,subprocess,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    files=set();base=ROOT/'results/uncertainty/development'
    for name in ['bispectrum-view-variance-v1','bispectrum-preferred-view-v1']:
        for ds in ['10028','10049']:assert json.loads((base/name/ds/'summary.json').read_text())['complete']
        files.update(p for p in (base/name).rglob('*') if p.is_file())
    name='bispectrum-view-risk-v1';assert json.loads((base/name/'summary.json').read_text())['complete']
    files.update(p for p in (base/name).rglob('*') if p.is_file())
    audit=ROOT/'research/uncertainty/reviews/view-variance-audit-01';manifest=json.loads((audit/'manifest.json').read_text())
    assert manifest['returncode']==0 and manifest['sources_unchanged'] and manifest['requested_model']=='claude-fable-5-1'
    files.update(p for p in audit.iterdir() if p.is_file())
    for prefix,names in {
        'src/fourier_splats':['uq_view_variance.py','uq_preferred_views.py','uq_view_risk.py','uq_moment_mc.py'],
        'tests':['test_view_variance.py','test_preferred_views.py','test_view_risk.py','test_moment_mc.py'],
        'scripts':['probe_uq_view_variance.py','probe_uq_preferred_views.py','probe_uq_view_risk.py',
            'verify_uq_view_variance.py','verify_uq_view_variance_sample.py','verify_uq_preferred_views.py','verify_uq_view_risk.py',
            'write_uq_view_variance_results.py','write_uq_view_followup_results.py','review_uq_view_variance.py',
            'package_uq_revision7.py','verify_uq_release_artifacts.py','build_paper.sh'],
        'research/uncertainty':['REPRODUCE-REVISION7.md','RELEASE-v0.7.4-dev.md','CURRENT-EVIDENCE.md'],
        'research/uncertainty/paired-power-v1':['VIEW-VARIANCE-PROTOCOL.md','VIEW-VARIANCE-RESULTS.md','VIEW-VARIANCE-PRIOR-ART.md',
            'PREFERRED-VIEW-PROTOCOL.md','PREFERRED-VIEW-RESULTS.md','VIEW-RISK-BASELINES-PROTOCOL.md','VIEW-RISK-BASELINES-RESULTS.md'],
        'provenance/uncertainty':['view-variance-independent-verification.json','view-variance-one-percent-verification.json',
            'preferred-views-independent-verification.json','view-risk-independent-verification.json','revision7-paper-verification.json',
            'view-variance-prior-art-downloads.json','view-risk-primary-downloads.json'],
        'paper':['main.tex','focused-main.tex','higher-moment-development.tex','view-variance-development.tex','uncertainty-references.bib']
    }.items():files.update(ROOT/prefix/name for name in names)
    for name in ['view-variance-initial-tests','bispectrum-view-variance-10028-v1','bispectrum-view-variance-10049-v1',
        'view-variance-independent-verification','view-variance-one-percent-verification','preferred-views-sampler-tests',
        'preferred-views-10028-v1','preferred-views-10049-v1','preferred-views-independent-verification',
        'view-risk-baseline-tests','view-risk-independent-verification','view-followup-final-tests','revision7-paper-build','revision7-render']:
        files.add(ROOT/'logs/uncertainty'/f'{name}.log')
    rows=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]
    path=ROOT/'output/artifacts/v0.7.4-dev-viewing-uncertainty.tar.gz';mp=path.with_name(path.name[:-7]+'-manifest.json')
    if path.exists() or mp.exists():raise ValueError('Preserve existing artifacts')
    with tarfile.open(path,'w:gz',compresslevel=1) as archive:
        for p in sorted(files):archive.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
    expected={r['path']:r for r in rows};seen=set()
    with tarfile.open(path,'r|gz') as archive:
        for member in archive:
            assert member.isfile() and member.name in expected and member.name not in seen
            with archive.extractfile(member) as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
            assert member.size==expected[member.name]['bytes'] and digest==expected[member.name]['sha256'];seen.add(member.name)
    assert seen==set(expected)
    record=dict(complete=True,stream_verification_complete=True,bytes=path.stat().st_size,sha256=sha(path),files=rows,
        source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope='Conditional known-simulator method development; not calibrated experimental uncertainty or a favorable full review.',
        dependencies=[dict(archive='v0.7.2-dev-method-gates.tar.gz',sha256='7fbb5540b0d1287716a457df2fbdcd6d6dcd50d3ba5b5cc08234a050599d533b',scope='Exact maps, geometry and original orbit arrays'),
            dict(archive='v0.7.3-dev-invariant-moments.tar.gz',sha256='001d6f46c6a655327efa4fa7050600e8ea56b4cbec7315739bcb8dbbbda743d3',scope='Frozen moment directions and independent earlier calibration')])
    mp.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k!='files'}),flush=True);print('members',len(rows),flush=True)


if __name__=='__main__':main()
