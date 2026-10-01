#!/usr/bin/env python3
"""Immutable evidence increment; explicit inclusion excludes private review data."""
import hashlib
import json
import subprocess
import tarfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def package(label,files):
    files=sorted(set(files));out=ROOT/f'output/artifacts/v0.7.8-dev-method-gates-{label}.tar.gz'
    manifest=out.with_name(out.name[:-7]+'-manifest.json')
    if out.exists() or manifest.exists():raise ValueError('Preserve immutable artifacts')
    assert all(p.is_file() for p in files)
    assert not any(p.name=='response.jsonl' for p in files)
    records=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in files]
    with tarfile.open(out,'w:gz',compresslevel=1) as tf:
        for p in files:tf.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
    expected={r['path']:r for r in records};seen=set()
    with tarfile.open(out,'r|gz') as tf:
        for member in tf:
            assert member.isfile() and member.name in expected and member.name not in seen
            with tf.extractfile(member) as stream:actual=hashlib.file_digest(stream,'sha256').hexdigest()
            assert actual==expected[member.name]['sha256'] and member.size==expected[member.name]['bytes']
            seen.add(member.name)
    assert set(expected)==seen and out.stat().st_size<2_000_000_000
    dependencies=[]
    for pattern in ['v0.7.7-dev-information-diagnostics-*-manifest.json','v0.7.6-dev-score-calibration-*-manifest.json']:
        for p in sorted((ROOT/'output/artifacts').glob(pattern)):
            j=json.loads(p.read_text())
            dependencies.append(dict(archive=p.name.replace('-manifest.json','.tar.gz'),sha256=j['sha256']))
    assert len(dependencies)==8
    j=dict(complete=True,stream_verification_complete=True,part=label,bytes=out.stat().st_size,sha256=sha(out),
           source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),files=records,
           dependencies=dependencies,scope='Classical integration repair, failed population-method screen, focused review; no new paper or acceptance verdict',
           manuscript='Unchanged v0.7.6-dev PDF is distributed in its earlier release')
    manifest.write_text(json.dumps(j,indent=2)+'\n');print(label,len(records),out.stat().st_size,j['sha256'],flush=True)


def main():
    for ds in ['10028','10049','10076']:
        files=[]
        for study in ['catalog-pose-integration-v1','population-materiality-v1']:
            d=BASE/study/ds
            assert json.loads((d/'summary.json').read_text())['complete']
            assert json.loads((d/'independent-check.json').read_text())['complete']
            files.extend(p for p in d.rglob('*') if p.is_file())
        package(ds,files)
    files=[]
    for study in ['catalog-pose-integration-report-v1','population-materiality-report-v1',
                  'adaptive-pose-autopsy-v1','catalog-pose-autopsy-v1','fixed-bank-separation-v1',
                  'acquisition-background-description-v1','background-source-group-correction-v1',
                  'cahra-inplane-convention-v1']:
        files.extend(p for p in (BASE/study).rglob('*') if p.is_file())
    files.extend(p for p in (BASE/'catalog-pose-integration-v1/logs').glob('*.log'))
    for pattern in ['CATALOG-POSE-*.md','POPULATION-MATERIALITY-*.md','POPULATION-SCREEN-*.md',
                    'AMPLITUDE-UNBIASED-*.md','PLANTED-TRUTH-*.md','FIXED-BANK-*.md',
                    'ACQUISITION-BACKGROUND-*.md','BACKGROUND-SOURCE-GROUP-*.md',
                    'CAHRA-POSE-*.md','cahra-pose-*.json','CRYOBIFE-SYNTHETIC-*.md',
                    'cryobife-supplement-*.json','CRYOTWIN-UNCERTAINTY-*.md','cryotwin*-source.json']:
        files.extend((ROOT/'research/uncertainty').glob(pattern))
    files.extend(ROOT/'research/uncertainty'/n for n in ['RELEASE-v0.7.8-dev.md','REPRODUCE-POST-REVIEW4-GATES.md'])
    review=ROOT/'research/uncertainty/reviews/post-round04-method-consultation'
    files.extend(review/n for n in ['critique.md','assistant-text.md','prompt.txt','manifest.json',
        'response.json','response.md','public-events.jsonl','publication-audit.json'])
    names=['probe_catalog_pose_integration.py','verify_catalog_pose_integration.py','write_catalog_pose_integration_report.py',
           'screen_population_materiality.py','verify_population_materiality.py','write_population_materiality_report.py',
           'autopsy_pose_integration.py','report_fixed_bank_separation.py','retain_fable_consultation.py',
           'describe_acquisition_background.py','verify_background_source_groups.py','audit_cahra_inplane_conventions.py',
           'package_post_review4_gates.py']
    files.extend(ROOT/'scripts'/n for n in names)
    files.extend(ROOT/'src/fourier_splats'/n for n in ['uq_pose_importance.py','uq_pose_catalog.py','uq_population_screen.py'])
    files.extend(ROOT/'tests'/n for n in ['test_pose_importance.py','test_pose_catalog.py','test_population_screen.py','test_planted_importance_identity.py'])
    package('common',files)


if __name__=='__main__':main()
