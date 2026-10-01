#!/usr/bin/env python3
"""Package immutable post-review-4 diagnostics, without a new paper claim."""
import hashlib
import json
import subprocess
import tarfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
DEPENDENCIES={
    '10028':'47b1860f3506f879169dce1cb04d1a281d811461298deca5b35b6fc4d160a59e',
    '10049':'c282a9b1f318c4ee00be4e5a8bfe96648c323272f64a1c0ca9f81a49bae05991',
    '10076':'1808b0b99a5165ddf77669126697aeaf13678bea8a45be7ebd8a5a9765c8da9a'}


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def package(label,files):
    files=sorted(set(files));path=ROOT/f'output/artifacts/v0.7.7-dev-information-diagnostics-{label}.tar.gz'
    manifest=path.with_name(path.name[:-7]+'-manifest.json')
    if path.exists() or manifest.exists():raise ValueError('Preserve immutable artifact attempts')
    rows=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in files]
    with tarfile.open(path,'w:gz',compresslevel=1) as archive:
        for p in files:archive.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
    expected={r['path']:r for r in rows};seen=set()
    with tarfile.open(path,'r|gz') as archive:
        for member in archive:
            assert member.isfile() and member.name in expected and member.name not in seen
            with archive.extractfile(member) as f:checksum=hashlib.file_digest(f,'sha256').hexdigest()
            assert member.size==expected[member.name]['bytes'] and checksum==expected[member.name]['sha256']
            seen.add(member.name)
    assert seen==set(expected) and path.stat().st_size<2_000_000_000
    result=dict(complete=True,stream_verification_complete=True,part=label,bytes=path.stat().st_size,sha256=sha(path),files=rows,
        source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope='Post-review-4 information/nuisance diagnostics and a failed numerical integration gate; no new manuscript or acceptance result.',
        dependencies=[dict(archive=f'v0.7.6-dev-score-calibration-{ds}.tar.gz',sha256=s) for ds,s in DEPENDENCIES.items()],
        manuscript='Unchanged v0.7.6-dev PDF, distributed in that earlier release; deliberately not duplicated here.')
    manifest.write_text(json.dumps(result,indent=2)+'\n')
    print(label,len(rows),path.stat().st_size,result['sha256'],flush=True)


def main():
    for ds in ['10028','10049','10076']:
        files=[]
        for study in ['matched-information-ledger-v1','matched-haar-information-v1','adaptive-pose-integration-v1']:
            folder=BASE/study/ds;j=json.loads((folder/'summary.json').read_text());assert j['complete']
            files.extend(p for p in folder.rglob('*') if p.is_file())
        assert json.loads((BASE/'adaptive-pose-integration-v1'/ds/'independent-check.json').read_text())['complete']
        files.append(BASE/'background-spectrum-inventory-v2'/(ds+'-moments.npz'))
        package(ds,files)
    common=[]
    for study in ['matched-information-verification-v1','matched-information-report-v1','matched-information-report-v2',
        'matched-event-decomposition-v1','adaptive-pose-integration-report-v1','background-spectrum-inventory-v2',
        'stack-nuisance-inventory-v1','population-author-replay-v1']:
        common.extend(p for p in (BASE/study).rglob('*') if p.is_file() and p.suffix!='.npz')
    for pattern in ['MATCHED-*.md','SAVED-EVENT-*.md','ADAPTIVE-POSE-*.md','BACKGROUND-SPECTRUM-*.md',
        'STACK-NUISANCE-*.md','POPULATION-BASELINE-*.md','POPULATION-UQ-*.md','CAHRA-V2-*.md',
        'POSE-INTEGRATION-*.md','ORIENTATION-PRIOR-*.md','NUMERICAL-LIKELIHOOD-*.md']:
        common.extend((ROOT/'research/uncertainty').glob(pattern))
    for name in ['REPRODUCE-POST-REVIEW4.md','RELEASE-v0.7.7-dev.md']:
        common.append(ROOT/'research/uncertainty'/name)
    for name in ['analyze_matched_information.py','analyze_matched_haar_likelihood.py','verify_matched_information.py',
        'write_matched_information_report.py','report_saved_event_decomposition.py','probe_adaptive_pose_integration.py',
        'verify_adaptive_pose_integration.py','write_adaptive_pose_integration_report.py',
        'inventory_background_spectra.py','write_background_spectrum_report.py','inventory_stack_nuisances.py',
        'verify_stack_nuisance_inventory.py','write_stack_nuisance_inventory.py','replay_population_baseline.py',
        'check_population_baseline_replay.py','package_post_review4_diagnostics.py']:
        common.append(ROOT/'scripts'/name)
    for name in ['uq_pose_importance.py','uq_pose_cone_probe.py','uq_continuous.py','uq_data.py','uq_bispectrum.py']:
        common.append(ROOT/'src/fourier_splats'/name)
    for name in ['test_pose_importance.py','test_matched_haar_information.py']:
        common.append(ROOT/'tests'/name)
    for pattern in ['matched-information-*.log','adaptive-pose-*.log','background-spectrum-*.log']:
        common.extend((ROOT/'logs/uncertainty').glob(pattern))
    package('common',common)


if __name__=='__main__':main()
