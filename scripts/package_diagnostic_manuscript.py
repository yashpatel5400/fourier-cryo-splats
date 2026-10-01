#!/usr/bin/env python3
"""Package the diagnostic manuscript using an explicit public-file allowlist."""
import hashlib
import json
import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results/uncertainty/development'


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    archive = ROOT / 'output/artifacts/v0.7.9-dev-diagnostic-manuscript.tar.gz'
    manifest = archive.with_name(archive.name[:-7] + '-manifest.json')
    if archive.exists() or manifest.exists():
        raise ValueError('Preserve existing release artifacts')
    files = {ROOT / p for p in [
        'README.md', 'paper/main.tex', 'paper/diagnostic-main.tex',
        'paper/diagnostic-appendix.tex', 'paper/diagnostic-references.bib',
        'paper/references.bib', 'paper/uncertainty-references.bib',
        'paper/survey-references.bib', 'paper/icml2026.sty', 'paper/icml2026.bst',
        'scripts/build_paper.sh', 'scripts/write_diagnostic_paper_tables.py',
        'scripts/package_diagnostic_manuscript.py',
        'output/pdf/fourier-cryo-splats.pdf',
        'provenance/uncertainty/diagnostic-manuscript-verification.json',
        'logs/uncertainty/diagnostic-manuscript-build.log',
        'research/uncertainty/REPRODUCE-DIAGNOSTIC-MANUSCRIPT.md',
        'research/uncertainty/RELEASE-v0.7.9-dev.md',
        'research/uncertainty/REPRODUCE-POST-REVIEW4-GATES.md',
        'research/uncertainty/CURRENT-EVIDENCE.md',
        'research/uncertainty/reviews/response-to-round-04-development.md',
        'results/uncertainty/development/catalog-pose-integration-report-v1/final-integration-gate.pdf',
    ]}
    files.update(ROOT / f'paper/tables/diagnostic-{name}.tex'
                 for name in ['separation', 'integration', 'population'])
    pages = sorted((ROOT / 'output/qa/diagnostic-rewrite').glob('final-page-*.png'))
    assert len(pages) == 16
    files.update(pages)
    for ds in ['10028', '10049', '10076']:
        for name in ['matched-information-ledger-v1', 'matched-haar-information-v1',
                     'adaptive-pose-integration-v1', 'catalog-pose-integration-v1',
                     'population-materiality-v1']:
            files.add(BASE / name / ds / 'summary.json')
        for name in ['catalog-pose-integration-v1', 'population-materiality-v1']:
            files.add(BASE / name / ds / 'independent-check.json')
    dependencies = []
    for path in sorted((ROOT / 'output/artifacts').glob('v0.7.8-dev-*-manifest.json')):
        record = json.loads(path.read_text())
        dependencies.append(dict(archive=path.name.replace('-manifest.json', '.tar.gz'),
                                 sha256=record['sha256'], bytes=record['bytes']))
        files.add(path)
    assert len(dependencies) == 4
    assert all(p.is_file() and not p.is_symlink() for p in files)
    rows = [dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, sha256=sha(p))
            for p in sorted(files)]
    with tarfile.open(archive, 'w:gz', compresslevel=6) as output:
        for path in sorted(files):
            output.add(path, arcname=str(path.relative_to(ROOT)), recursive=False)
    expected = {r['path']: r for r in rows}
    seen = set()
    with tarfile.open(archive, 'r|gz') as stream:
        for member in stream:
            assert member.isfile() and member.name in expected and member.name not in seen
            with stream.extractfile(member) as data:
                digest = hashlib.file_digest(data, 'sha256').hexdigest()
            assert digest == expected[member.name]['sha256']
            assert member.size == expected[member.name]['bytes']
            seen.add(member.name)
    assert seen == set(expected)
    result = dict(complete=True, stream_verification_complete=True, files=rows,
                  archive=archive.name, bytes=archive.stat().st_size, sha256=sha(archive),
                  source_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                                      cwd=ROOT, text=True).strip(),
                  dependencies=dependencies,
                  scope='Diagnostic manuscript; no new experimental calibration or acceptance verdict.')
    manifest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'files'}, indent=2))
    print('Verified members:', len(rows))


if __name__ == '__main__':
    main()
