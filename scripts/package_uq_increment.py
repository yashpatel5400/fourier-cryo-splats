#!/usr/bin/env python3
"""Package an explicit, hash-locked research artifact specification.

This creates local files only. It never uploads, deletes, extracts an archive,
or chooses experiment results. The spec lists every member and its completed
record owner, with expected SHA-256 digests for both.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def local(path):
    resolved = (ROOT / path).resolve(strict=True)
    resolved.relative_to(ROOT)
    if not resolved.is_file() or Path(path).is_absolute() or '..' in Path(path).parts:
        raise ValueError(f'Not a repository-relative file: {path}')
    return resolved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('spec', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    spec = json.loads(args.spec.read_text())
    archive = args.output.resolve()
    if not archive.name.endswith('.tar.gz'):
        raise ValueError('Use a .tar.gz output path')
    manifest = archive.with_name(archive.name[:-7] + '-manifest.json')
    if archive.exists() or manifest.exists():
        raise FileExistsError('Preserve existing artifacts; choose a new output')
    archive.parent.mkdir(parents=True, exist_ok=True)
    entries = sorted(spec['files'], key=lambda row: row['path'])
    if not entries or len({row['path'] for row in entries}) != len(entries):
        raise ValueError('Empty or duplicate member list')
    owners = {}
    members = []
    for entry in entries:
        path = local(entry['path'])
        if digest(path) != entry['sha256']:
            raise ValueError(f"Member hash changed: {entry['path']}")
        owner = local(entry['owner'])
        owner_hash = digest(owner)
        if owner_hash != entry['owner_sha256']:
            raise ValueError(f"Owner hash changed: {entry['owner']}")
        record = json.loads(owner.read_text())
        if not record.get('complete') or record.get('error') or record.get('scientific_run_complete') is False:
            raise ValueError(f"Owner is not a completed successful record: {entry['owner']}")
        owners[entry['owner']] = owner_hash
        members.append(dict(path=entry['path'], sha256=entry['sha256'],
                            bytes=path.stat().st_size, owner=entry['owner']))
    # Fixed metadata makes identical member bytes produce identical archives.
    # Exclusive creation preserves both successful and partial historical outputs.
    with archive.open('xb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode='w|', format=tarfile.PAX_FORMAT) as tar:
                for member in members:
                    path = local(member['path'])
                    info = tarfile.TarInfo(member['path'])
                    info.size = member['bytes']
                    info.mode = 0o644
                    with path.open('rb') as handle:
                        tar.addfile(info, handle)
    expected = {row['path']: row for row in members}
    seen = set()
    with tarfile.open(archive, mode='r|gz') as tar:
        for member in tar:
            if not member.isfile() or member.name not in expected or member.name in seen:
                raise ValueError(f'Unexpected archive member: {member.name}')
            with tar.extractfile(member) as handle:
                actual = hashlib.file_digest(handle, 'sha256').hexdigest()
            if actual != expected[member.name]['sha256'] or member.size != expected[member.name]['bytes']:
                raise ValueError(f'Archive verification failed: {member.name}')
            seen.add(member.name)
    if seen != set(expected):
        raise ValueError('Archive member set differs from specification')
    # Detect a concurrent owner change before declaring the snapshot complete.
    for owner, expected_hash in owners.items():
        if digest(local(owner)) != expected_hash:
            raise ValueError(f'Owner changed during packaging: {owner}')
    result = dict(complete=True, version=spec['version'], archive=archive.name,
        sha256=digest(archive), bytes=archive.stat().st_size, files=members,
        completed_record_owners=owners, spec_sha256=digest(args.spec),
        packaging_script_sha256=digest(Path(__file__)),
        source_git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        required_prior_arrays=spec['required_prior_arrays'], scope=spec['scope'],
        reproduction=spec['reproduction'], stream_verification_complete=True)
    with manifest.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(archive=str(archive), manifest=str(manifest),
                         sha256=result['sha256'], bytes=result['bytes'], members=len(members))))


if __name__ == '__main__':
    main()
