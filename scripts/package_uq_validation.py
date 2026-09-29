#!/usr/bin/env python3
"""Package exact validation weights/poses; original trained models remain v0.2."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', default='v0.3.0-dev')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9.-]+', args.version):
        raise ValueError('Simple version identifier required')
    for name, count in [('continuous-v1', 96), ('continuous-moments-v2', 48)]:
        summary = json.loads((ROOT/f'results/uncertainty/confirmation/{name}/summary/summary.json').read_text())
        if summary['audit_settings'] != count or not summary['status'].startswith('complete'):
            raise AssertionError('Cannot package an incomplete frozen study')
    folders = [ROOT/'results/uncertainty/development'/name for name in [
        'continuous-optimized', 'continuous-quadrature-optimized',
        'continuous-quadrature-1024particles', 'continuous-quadrature-high-band-probe',
        'continuous-quadrature-high-band-additional', 'continuous-pose-adversaries']]
    folders += [ROOT/'results/uncertainty/confirmation'/name for name in [
        'continuous-v1', 'continuous-moments-v2', 'prediction-v1']]
    files = sorted({p for folder in folders for p in folder.rglob('*.npz')})
    if not files:
        raise AssertionError('No validation arrays found')
    out = ROOT/'output/artifacts'; out.mkdir(parents=True, exist_ok=True)
    archive = out/f'uncertainty-validation-arrays-{args.version}.tar.gz'
    manifest_path = ROOT/'provenance/uncertainty'/f'validation-arrays-{args.version}.json'
    if archive.exists() or manifest_path.exists():
        raise RuntimeError('Do not overwrite published or prepared artifacts')
    entries = [{'path': str(p.relative_to(ROOT)), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in files]
    with archive.open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as gz, tarfile.open(fileobj=gz, mode='w') as tar:
        for path in files:
            info = tar.gettarinfo(str(path), arcname=str(path.relative_to(ROOT)))
            info.mtime = 0; info.uid = 0; info.gid = 0; info.uname = ''; info.gname = ''
            with path.open('rb') as content:
                tar.addfile(info, content)
    manifest = {'version': args.version, 'archive': archive.name, 'bytes': archive.stat().st_size,
                'sha256': sha(archive), 'files': entries,
                'source_git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'scope': 'Exact selected estimator weights, nonlinear stress poses and frozen prediction errors; no particle pixels or third-party PDFs.',
                'required_prior_checkpoint': 'https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.2.0-dev',
                'reproduction': 'Extract at repository root; fetch original particle selections/maps using recorded protocols.'}
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({k: manifest[k] for k in ['archive', 'bytes', 'sha256']}, indent=2))


if __name__ == '__main__':
    main()
