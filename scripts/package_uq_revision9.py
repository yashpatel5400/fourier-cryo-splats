#!/usr/bin/env python3
"""Split complete score/allocation evidence into sub-2GB per-stack bundles."""
import hashlib,json,subprocess,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def package(label,files):
    files=sorted(files);path=ROOT/f'output/artifacts/v0.7.6-dev-score-calibration{label}.tar.gz';mp=path.with_name(path.name[:-7]+'-manifest.json')
    if path.exists() or mp.exists():raise ValueError('Preserve existing artifacts')
    rows=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in files]
    with tarfile.open(path,'w:gz',compresslevel=1) as archive:
        for p in files:archive.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
    expected={r['path']:r for r in rows};seen=set()
    with tarfile.open(path,'r|gz') as archive:
        for m in archive:
            assert m.isfile() and m.name in expected and m.name not in seen
            with archive.extractfile(m) as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
            assert m.size==expected[m.name]['bytes'] and digest==expected[m.name]['sha256'];seen.add(m.name)
    assert seen==set(expected) and path.stat().st_size<2_000_000_000
    record=dict(complete=True,stream_verification_complete=True,bytes=path.stat().st_size,sha256=sha(path),files=rows,
        source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope='Known-simulator classical score and allocation diagnostics; not experimentally calibrated local density coverage.',
        part=label or 'common',required_parts=['common','-10028','-10049','-10076'],
        dependencies=[dict(archive='v0.7.5-dev-candidate-scores.tar.gz',sha256='32170df45f4aa2833a6e5cfd541881431034450ce5a7d8fc4af1f803a7ab002b',
            scope='Actual fitted Gaussian maps and preceding development; earlier dependencies as recorded in that manifest.')])
    mp.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k!='files'}),flush=True);print('members',len(rows),flush=True)
def main():
    names=['candidate-fisher-score-v1','candidate-replica-allocation-v1'];common=set()
    for ds in ['10028','10049','10076']:
        files=set()
        for name in names:
            sp=BASE/name/ds/'summary.json';r=json.loads(sp.read_text());assert r['complete'];assert sha(ROOT/r['array']['path'])==r['array']['sha256']
            files.update(p for p in (BASE/name/ds).rglob('*') if p.is_file())
        package('-'+ds,files)
    for name in names:common.update(p for p in (BASE/name).iterdir() if p.is_file())
    for prefix,filenames in {
        'src/fourier_splats':['uq_fisher_score.py','uq_candidate_score.py','uq_view_variance.py','uq_view_risk.py','uq_moment_mc.py','uq_bispectrum.py','uq_continuous.py'],
        'tests':['test_fisher_score.py','test_candidate_score.py','test_view_variance.py','test_view_risk.py','test_moment_mc.py','test_preferred_views.py'],
        'scripts':['probe_uq_fisher_scores.py','probe_uq_replica_allocation.py','verify_uq_fisher_scores.py','verify_uq_fisher_comparisons.py',
            'verify_uq_replica_allocation.py','write_uq_fisher_score_results.py','write_uq_replica_allocation_results.py',
            'package_uq_revision9.py','verify_uq_release_artifacts.py','build_paper.sh'],
        'research/uncertainty':['REPRODUCE-REVISION9.md','RELEASE-v0.7.6-dev.md','CURRENT-EVIDENCE.md'],
        'research/uncertainty/paired-power-v1':['FISHER-SCORE-PROTOCOL.md','FISHER-SCORE-COMPARISON-INTERVALS.md','FISHER-SCORE-DERIVATION.md',
            'FISHER-SCORE-RESULTS.md','NESTED-RISK-READING.md','REPLICA-ALLOCATION-PROTOCOL.md','REPLICA-ALLOCATION-RESULTS.md'],
        'provenance/uncertainty':['fisher-score-independent-verification.json','fisher-comparison-independent-verification.json',
            'replica-allocation-independent-verification.json','nested-risk-reading-downloads.json','revision9-paper-verification.json'],
        'paper':['candidate-score-development.tex','fisher-score-development.tex','view-variance-development.tex','main.tex','focused-main.tex','uncertainty-references.bib'],
        'paper/figures':['fisher-score-projections.pdf','fisher-score-projections.png'],
        'logs/uncertainty':['fisher-score-initial-tests.log','fisher-followup-tests.log','fisher-scores-10028-v1.log','fisher-scores-10049-v1.log','fisher-scores-10076-v1.log',
            'fisher-score-independent-verification.log','fisher-comparison-independent-verification.log','fisher-score-report.log',
            'replica-allocation-10028-v1.log','replica-allocation-10049-v1.log','replica-allocation-10076-v1.log',
            'replica-allocation-independent-verification.log','replica-allocation-report.log','revision9-paper-build.log','revision9-render.log']
    }.items():common.update(ROOT/prefix/name for name in filenames)
    package('',common)
if __name__=='__main__':main()
