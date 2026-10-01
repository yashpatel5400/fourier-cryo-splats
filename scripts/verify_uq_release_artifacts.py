#!/usr/bin/env python3
"""Read-only verification of a published archive and every declared member."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile


def verify(archive,manifest):
    record=json.loads(manifest.read_text())
    if not record.get('complete') or not record.get('stream_verification_complete'):
        raise ValueError('Manifest does not describe a completed verified package')
    with archive.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
    if archive.stat().st_size!=record['bytes'] or actual!=record['sha256']:
        raise ValueError('Archive size/hash differs from manifest')
    expected={r['path']:r for r in record['files']}
    if len(expected)!=len(record['files']):raise ValueError('Duplicate manifest path')
    seen=set()
    with tarfile.open(archive,mode='r|gz') as stream:
        for member in stream:
            if (not member.isfile() or member.name not in expected or member.name in seen
                or Path(member.name).is_absolute() or '..' in Path(member.name).parts):
                raise ValueError('Unexpected or unsafe archive member: '+member.name)
            row=expected[member.name]
            with stream.extractfile(member) as f:h=hashlib.file_digest(f,'sha256').hexdigest()
            if member.size!=row['bytes'] or h!=row['sha256']:
                raise ValueError('Member size/hash differs: '+member.name)
            seen.add(member.name)
    if seen!=set(expected):raise ValueError('Missing members')
    return dict(archive=str(archive),sha256=actual,verified_members=len(seen),bytes=record['bytes'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archives',nargs='+',type=Path)
    args=parser.parse_args()
    for path in args.archives:
        if not path.name.endswith('.tar.gz'):raise ValueError('Expected .tar.gz archive')
        manifest=path.with_name(path.name[:-7]+'-manifest.json')
        print(json.dumps(verify(path,manifest)),flush=True)


if __name__=='__main__':main()
