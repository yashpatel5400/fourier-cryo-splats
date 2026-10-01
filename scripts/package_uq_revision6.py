#!/usr/bin/env python3
"""Package the complete invariant-moment development increment."""
import hashlib,json,subprocess,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    files=set();base=ROOT/'results/uncertainty/development'
    for name in ['bispectrum-hull-v1','bispectrum-adversarial-pose-v1',
                 'bispectrum-global-bound-v1','bispectrum-monte-carlo-v1']:
        summary=json.loads((base/name/'summary.json').read_text())
        assert summary['complete']
        files.update(p for p in (base/name).rglob('*') if p.is_file())
    for prefix,names in {
        'src/fourier_splats':['uq_bispectrum.py','uq_moment_rotation_bound.py','uq_moment_mc.py'],
        'tests':['test_bispectrum.py','test_moment_rotation_bound.py','test_moment_mc.py'],
        'scripts':['probe_uq_bispectrum.py','probe_uq_bispectrum_poses.py','probe_uq_moment_global_bound.py',
            'probe_uq_moment_mc.py','verify_uq_bispectrum.py','verify_uq_moment_mc.py',
            'write_uq_bispectrum_results.py','write_uq_moment_mc_results.py',
            'package_uq_revision6.py','verify_uq_release_artifacts.py','build_paper.sh'],
        'research/uncertainty':['ALIGNMENT-LITERATURE-REVIEW3.md','alignment-review3-source-reading.json',
            'REPRODUCE-REVISION6.md','CURRENT-EVIDENCE.md','SURVEY.md','reading-list.tsv'],
        'research/uncertainty/paired-power-v1':['BISPECTRUM-PROTOCOL.md','BISPECTRUM-POSE-PROTOCOL.md',
            'BISPECTRUM-RESULTS.md','BISPECTRUM-GLOBAL-BOUND-PROTOCOL.md','BISPECTRUM-GLOBAL-BOUND-RESULTS.md',
            'MONTE-CARLO-VIEW-LAW-PROTOCOL.md','MONTE-CARLO-VIEW-LAW-RESULTS.md'],
        'provenance/uncertainty':['alignment-moments-followup-downloads.json','bispectrum-independent-verification.json',
            'moment-mc-independent-verification.json','revision6-paper-verification.json','revision6-mc-paper-verification.json'],
        'paper':['main.tex','focused-main.tex','focused-survey.tex','higher-moment-development.tex',
            'uncertainty-references.bib','tables/moment-mc-development.tex']
    }.items():files.update(ROOT/prefix/name for name in names)
    logs=['bispectrum-initial-tests','bispectrum-prefreeze-tests','bispectrum-hull-v1',
        'bispectrum-adversarial-pose-v1','bispectrum-independent-verification','moment-rotation-bound-tests',
        'bispectrum-global-bound-v1','moment-mc-initial-tests','bispectrum-monte-carlo-v1',
        'moment-mc-independent-verification','revision6-paper-build','revision6-paper-build-final',
        'revision6-render','revision6-mc-paper-build','revision6-mc-render']
    files.update(ROOT/'logs/uncertainty'/f'{name}.log' for name in logs)
    rows=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]
    path=ROOT/'output/artifacts/v0.7.3-dev-invariant-moments.tar.gz'
    manifest=path.with_name(path.name[:-7]+'-manifest.json')
    if path.exists() or manifest.exists():raise ValueError('Preserve existing release artifacts')
    with tarfile.open(path,'w:gz',compresslevel=1) as archive:
        for p in sorted(files):archive.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
    expected={r['path']:r for r in rows};seen=set()
    with tarfile.open(path,'r|gz') as archive:
        for member in archive:
            assert member.isfile() and member.name in expected and member.name not in seen
            with archive.extractfile(member) as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
            assert member.size==expected[member.name]['bytes'] and digest==expected[member.name]['sha256']
            seen.add(member.name)
    assert seen==set(expected)
    record=dict(complete=True,stream_verification_complete=True,bytes=path.stat().st_size,sha256=sha(path),
        files=rows,source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope='Post-review method development under explicit oracle/known-noise/bounded-view assumptions; not experimental density coverage or acceptance.',
        dependencies=['Tagged source and pinned Python environment.',
            'v0.7.2-dev-method-gates.tar.gz (sha256 7fbb5540b0d1287716a457df2fbdcd6d6dcd50d3ba5b5cc08234a050599d533b) supplies exact orbit arrays, three EMDB maps and preflight geometry.'])
    manifest.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='files'}),flush=True);print('members',len(rows),flush=True)


if __name__=='__main__':main()
