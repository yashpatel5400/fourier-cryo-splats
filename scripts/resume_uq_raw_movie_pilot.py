#!/usr/bin/env python3
"""Resume the interrupted, preselected movie without replacing failed evidence."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import requests
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]


def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    script=Path(__file__)
    if subprocess.check_output(['git','show',f'HEAD:{script.relative_to(ROOT)}'],cwd=ROOT)!=script.read_bytes():
        raise ValueError('Commit the transfer recovery before execution')
    canonical=ROOT/'provenance/uncertainty/raw-movie-pilot-v1.json'
    old=json.loads(canonical.read_text())
    if old['complete'] or not old.get('error'):
        raise ValueError('Only the recorded incomplete transfer may be resumed')
    directory=ROOT/'data/uncertainty/raw-movie-pilot-v1'
    prefix=directory/'004_movie.mrcs.partial'
    tail=directory/'004_movie.mrcs.resume-tail'
    final=directory/'004_movie.mrcs'
    archive=canonical.with_name('raw-movie-pilot-download-attempt-01.json')
    record_path=canonical.with_name('raw-movie-pilot-resume-v1.json')
    if any(p.exists() for p in [archive,tail,final,record_path]):
        raise RuntimeError('Preserve previous recovery outputs; do not overwrite')
    start=prefix.stat().st_size; total=old['expected_bytes']
    if start!=old['bytes_downloaded'] or not 1024<start<total:
        raise ValueError('Partial file size disagrees with saved failure')
    header=ROOT/'data/uncertainty/movie-feasibility/004-header.mrc'
    with prefix.open('rb') as f:
        if f.read(1024)!=header.read_bytes():raise ValueError('Partial movie header changed')
    archive.write_bytes(canonical.read_bytes())
    prefix_hash=sha(prefix);tick=time.perf_counter()
    record=dict(complete=False,offset=start,total_bytes=total,tail_bytes=0,
        previous_attempt=str(archive.relative_to(ROOT)),previous_attempt_sha256=sha(archive),
        prefix_sha256=prefix_hash,source_snapshot=source_snapshot(ROOT,script,[]),
        scope='Transport recovery only; the source movie and frozen diagnostic analysis are unchanged.')
    def save():record_path.write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        etag=old['response_headers']['ETag']
        with requests.get(old['url'],headers={'Range':f'bytes={start}-{total-1}',
                'If-Range':etag,'Accept-Encoding':'identity'},stream=True,timeout=(20,60)) as response:
            response.raise_for_status()
            record.update(status=response.status_code,final_url=response.url,
                response_headers={k:response.headers.get(k) for k in ['Content-Range','Content-Length','ETag','Last-Modified']})
            save()
            if (response.status_code!=206 or response.headers.get('Content-Range')!=f'bytes {start}-{total-1}/{total}'
                or int(response.headers.get('Content-Length',-1))!=total-start
                or response.headers.get('ETag')!=etag
                or response.headers.get('Last-Modified')!=old['response_headers']['Last-Modified']):
                raise ValueError('Remote range or source identity changed')
            with tail.open('xb') as f:
                for block in response.iter_content(4*1024*1024):
                    f.write(block);record['tail_bytes']+=len(block)
                    if record['tail_bytes']>total-start:raise ValueError('Oversized range')
                    save()
        if tail.stat().st_size!=total-start or sha(prefix)!=prefix_hash:
            raise ValueError('Transfer truncated or preserved prefix changed')
        with final.open('xb') as f:
            for p in [prefix,tail]:
                with p.open('rb') as source:shutil.copyfileobj(source,f,8*1024*1024)
        if final.stat().st_size!=total:raise ValueError('Assembled size differs')
        with final.open('rb') as f:
            if f.read(1024)!=header.read_bytes():raise ValueError('Assembled header differs')
        record.update(complete=True,tail_sha256=sha(tail),assembled_sha256=sha(final),seconds=time.perf_counter()-tick)
        save()
        updated=dict(old)
        updated.pop('error')
        updated.update(complete=True,path=str(final.relative_to(ROOT)),sha256=record['assembled_sha256'],
            bytes_downloaded=total,seconds=old['seconds']+record['seconds'],
            previous_failed_attempt=str(archive.relative_to(ROOT)),previous_failed_attempt_sha256=sha(archive),
            resume_record=str(record_path.relative_to(ROOT)),resume_record_sha256=sha(record_path),
            transfer_resumed=True)
        canonical.write_text(json.dumps(updated,indent=2)+'\n')
        print('Verified resumed movie',total,record['assembled_sha256'],flush=True)
    except Exception as exc:
        record.update(error=repr(exc),seconds=time.perf_counter()-tick);save();raise


if __name__=='__main__':main()
