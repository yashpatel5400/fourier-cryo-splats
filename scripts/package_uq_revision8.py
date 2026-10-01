#!/usr/bin/env python3
"""Package candidate-derived simulations and their actual Gaussian MRC inputs."""
import hashlib,json,subprocess,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    base=ROOT/'results/uncertainty/development/candidate-moment-score-v1'
    files={p for p in base.rglob('*') if p.is_file()}
    for ds in ['10028','10049','10076']:
        r=json.loads((base/ds/'summary.json').read_text());assert r['complete']
        p=ROOT/f'results/final/{ds}/gaussian-mean.mrc';assert sha(p)==r['input_hashes'][str(p.relative_to(ROOT))];files.add(p)
    assert json.loads((base/'classical-comparators.json').read_text())['complete']
    for prefix,names in {
        'src/fourier_splats':['uq_candidate_score.py','uq_view_variance.py','uq_view_risk.py','uq_moment_mc.py','uq_bispectrum.py','uq_continuous.py'],
        'tests':['test_candidate_score.py'],
        'scripts':['probe_uq_candidate_scores.py','analyze_uq_candidate_scores.py','verify_uq_candidate_scores.py',
            'write_uq_candidate_score_results.py','package_uq_revision8.py','verify_uq_release_artifacts.py','build_paper.sh'],
        'research/uncertainty':['REPRODUCE-REVISION8.md','RELEASE-v0.7.5-dev.md','CURRENT-EVIDENCE.md'],
        'research/uncertainty/paired-power-v1':['CANDIDATE-SCORE-PROTOCOL.md','CANDIDATE-SCORE-COMPARATORS.md','CANDIDATE-SCORE-RESULTS.md'],
        'provenance/uncertainty':['candidate-score-independent-verification.json','revision8-paper-verification.json'],
        'paper':['candidate-score-development.tex','view-variance-development.tex','main.tex','focused-main.tex','uncertainty-references.bib'],
        'paper/figures':['candidate-score-projections.pdf','candidate-score-projections.png'],
        'logs/uncertainty':['candidate-score-initial-tests.log','candidate-scores-10028-v1.log','candidate-scores-10049-v1.log','candidate-scores-10076-v1.log',
            'candidate-score-classical-comparators.log','candidate-score-independent-verification.log','revision8-paper-build.log','revision8-render.log']
    }.items():files.update(ROOT/prefix/name for name in names)
    rows=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]
    path=ROOT/'output/artifacts/v0.7.5-dev-candidate-scores.tar.gz';mp=path.with_name(path.name[:-7]+'-manifest.json')
    if path.exists() or mp.exists():raise ValueError('Preserve existing artifacts')
    with tarfile.open(path,'w:gz',compresslevel=1) as archive:
        for p in sorted(files):archive.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
    expected={r['path']:r for r in rows};seen=set()
    with tarfile.open(path,'r|gz') as archive:
        for m in archive:
            assert m.isfile() and m.name in expected and m.name not in seen
            with archive.extractfile(m) as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
            assert m.size==expected[m.name]['bytes'] and digest==expected[m.name]['sha256'];seen.add(m.name)
    assert seen==set(expected)
    record=dict(complete=True,stream_verification_complete=True,bytes=path.stat().st_size,sha256=sha(path),files=rows,
        source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope='Candidate-driven known-simulator tests; not regional occupancy coverage or calibrated experimental uncertainty.',
        dependencies=[dict(archive='v0.7.2-dev-method-gates.tar.gz',sha256='7fbb5540b0d1287716a457df2fbdcd6d6dcd50d3ba5b5cc08234a050599d533b'),
            dict(archive='v0.7.3-dev-invariant-moments.tar.gz',sha256='001d6f46c6a655327efa4fa7050600e8ea56b4cbec7315739bcb8dbbbda743d3'),
            dict(archive='v0.7.4-dev-viewing-uncertainty.tar.gz',sha256='1638defe3b687a0ef28482cacc2857144241e93c64cedc63e38d5eda84836ca4')])
    mp.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k!='files'}));print('members',len(rows))
if __name__=='__main__':main()
