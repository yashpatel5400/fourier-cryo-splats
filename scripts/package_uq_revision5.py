#!/usr/bin/env python3
"""Package complete method-gate arrays, dependencies and verified paper changes."""
import hashlib,json,subprocess,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    files=set();base=ROOT/'results/uncertainty/development'
    directories=['paired-power-enlarged-cone-v1','paired-covariance-finite-view-v1',
        'paired-covariance-full-frequency-v1','paired-covariance-shifted-v1',
        'paired-covariance-adversarial-10049-v1','paired-statistics-bound-corrections-v1',
        'paired-method-gate-summary-v1','folded-ridge-review3-v1','mixture-validation-preflight-v2']
    directories += [f'paired-power-finite-view-v{i}' for i in [1,2,3,4]]
    for directory in directories:
        files.update(p for p in (base/directory).rglob('*') if p.is_file())
    for ds,emd in [('10028','2660'),('10049','6487'),('10076','8434')]:
        files.add(ROOT/f'data/uncertainty/references/emd_{emd}.map')
        files.add(base/f'local-alignment-calibration-v1/{ds}/generators.npz')
        files.add(base/f'refitting-bias-reanalysis-v1/{ds}/true-pose-fits.json')
    for directory in ['research/uncertainty/paired-power-v1','research/uncertainty/reviews/paired-statistics-audit-01']:
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    # Source lives in the Git tag as well; include every module that implements
    # the new computations and every dedicated test/runner/report generator.
    for prefix,names in {
        'src/fourier_splats':['uq_paired_power.py','uq_power_cone.py','uq_power_cone_witness.py',
            'uq_paired_covariance.py','uq_covariance_cone_bounds.py','uq_folded_ridge.py','uq_shift_moments.py','uq_pose_cone_probe.py'],
        'tests':['test_paired_power.py','test_power_cone.py','test_power_cone_witness.py','test_paired_covariance.py',
            'test_covariance_cone_bounds.py','test_folded_ridge.py','test_shift_moments.py','test_pose_cone_probe.py'],
        'scripts':['probe_uq_enlarged_power_cone.py','verify_uq_power_cone_witnesses.py','probe_uq_paired_covariance.py',
            'review_uq_paired_statistics.py','probe_uq_full_covariance.py','probe_uq_shifted_covariance.py',
            'probe_uq_adversarial_pose.py','probe_uq_folded_ridge.py','correct_uq_paired_bound_diagnostics.py',
            'verify_uq_revision5_diagnostics.py','verify_uq_adversarial_witnesses.py','write_uq_folded_ridge_results.py',
            'write_uq_method_gate_results.py','package_uq_revision5.py','verify_uq_release_artifacts.py','build_paper.sh'],
        'research/uncertainty':['FOLDED-RIDGE-REVIEW3-PROTOCOL.md','FOLDED-RIDGE-REVIEW3-RESULTS.md','REPRODUCE-REVISION5.md'],
        'provenance/uncertainty':['enlarged-power-cone-verification.json','revision5-diagnostic-verification.json',
            'adversarial-pose-direct-verification.json','revision5-paper-verification.json'],
        'paper':['main.tex','focused-main.tex','focused-theory.tex','focused-experiments.tex','tables/folded-ridge-review3.tex']
    }.items():
        files.update(ROOT/prefix/name for name in names)
    log_names=['paired-covariance-primitive-tests','paired-statistics-audit-01','folded-ridge-tests-initial-failures',
        'folded-ridge-tests','folded-ridge-review3-v1','covariance-cone-bounds-tests','paired-statistics-post-audit-tests',
        'paired-statistics-bound-corrections-v1','shift-moments-tests','pose-cone-probe-tests',
        'paired-covariance-finite-view-v1','paired-covariance-full-frequency-v1','paired-covariance-shifted-v1',
        'paired-covariance-adversarial-10049-v1','revision5-diagnostic-verification','adversarial-pose-direct-verification',
        'revision5-focused-tests','revision5-paper-build','revision5-render']
    files.update(ROOT/'logs/uncertainty'/f'{name}.log' for name in log_names)
    rows=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]
    path=ROOT/'output/artifacts/v0.7.2-dev-method-gates.tar.gz'
    manifest=path.with_name(path.name[:-7]+'-manifest.json')
    if path.exists() or manifest.exists():raise ValueError('Preserve prior release package')
    with tarfile.open(path,'w:gz',compresslevel=1) as archive:
        for p in sorted(files):archive.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
    expected={r['path']:r for r in rows};seen=set()
    with tarfile.open(path,'r|gz') as archive:
        for member in archive:
            assert member.isfile() and member.name in expected and member.name not in seen
            with archive.extractfile(member) as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
            assert digest==expected[member.name]['sha256'] and member.size==expected[member.name]['bytes'];seen.add(member.name)
    assert seen==set(expected)
    record=dict(complete=True,stream_verification_complete=True,bytes=path.stat().st_size,sha256=sha(path),files=rows,
        source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope='Post-review development diagnostics and falsification outcomes, not a validated experimental uncertainty method.',
        dependencies=['Tagged source and pinned Python environment.',
            'Saved-array revision-5 replays use this bundle. Regenerating the older initial orbit/preflight from particle metadata requires the v0.7.0 dependencies and source inputs named in its original record.'])
    manifest.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='files'}),flush=True);print('members',len(rows),flush=True)


if __name__=='__main__':main()
