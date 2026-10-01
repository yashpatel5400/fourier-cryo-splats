#!/usr/bin/env python3
"""Download one declared raw movie with streaming size/header/provenance checks."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time
import requests
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT/'research/uncertainty/RAW-MOVIE-PILOT-PROTOCOL.md'
URL = 'https://ftp.ebi.ac.uk/empiar/world_availability/10028/data/Micrographs/Micrographs_part2/004_movie.mrcs'
EXPECTED_BYTES = 1073742848


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    for p in [Path(__file__), PROTOCOL]:
        if subprocess.check_output(['git', 'show', f'HEAD:{p.relative_to(ROOT)}'], cwd=ROOT) != p.read_bytes():
            raise ValueError('Commit selected movie/procedure before pixels')
    out = ROOT/'data/uncertainty/raw-movie-pilot-v1'; out.mkdir(parents=True, exist_ok=True)
    path = out/'004_movie.mrcs'; partial = path.with_suffix('.mrcs.partial')
    report_path = ROOT/'provenance/uncertainty/raw-movie-pilot-v1.json'
    if path.exists() or partial.exists() or report_path.exists():
        raise RuntimeError('Preserve prior movie download/partial/provenance')
    inputs = [ROOT/'research/uncertainty/splits/10028.csv',
        ROOT/'research/uncertainty/confirmation/prediction-v1/10028-selection.csv',
        ROOT/'research/uncertainty/confirmation/noise-calibration-v1/10028-selection.csv']
    used = set()
    for p in inputs:
        used.update(r['source_group'] for r in csv.DictReader(p.open()))
    if 'MRC_1901/004_movie.mrcs' in used:
        raise ValueError('Selected group was used')
    header = ROOT/'data/uncertainty/movie-feasibility/004-header.mrc'
    record = dict(complete=False, url=URL, expected_bytes=EXPECTED_BYTES, bytes_downloaded=0,
        source_snapshot=source_snapshot(ROOT, Path(__file__), [str(PROTOCOL.relative_to(ROOT))]),
        input_hashes={str(p.relative_to(ROOT)): sha(p) for p in inputs+[header]},
        excluded_groups=len(used), proposed_source_group='MRC_1901/004_movie.mrcs',
        identity_scope='Group mapping supported by counts/date/nearby recentered coordinates, not an exact identity-verified confirmatory cohort.')
    def save():
        report_path.write_text(json.dumps(record, indent=2)+'\n')
    save(); start = time.perf_counter()
    try:
        with requests.get(URL, stream=True, timeout=(20, 60)) as response:
            response.raise_for_status()
            record['response_headers'] = {k: response.headers.get(k) for k in ['Content-Length', 'Content-Type', 'ETag', 'Last-Modified']}
            record['final_url'] = response.url; save()
            if int(response.headers.get('Content-Length', -1)) != EXPECTED_BYTES:
                raise ValueError('Archive file length changed')
            h = hashlib.sha256()
            with partial.open('wb') as f:
                for block in response.iter_content(8*1024*1024):
                    f.write(block); h.update(block); record['bytes_downloaded'] += len(block)
                    if record['bytes_downloaded'] > EXPECTED_BYTES:
                        raise ValueError('Oversized movie')
                    if record['bytes_downloaded'] % (128*1024*1024) == 0:
                        save(); print('bytes', record['bytes_downloaded'], 'seconds', time.perf_counter()-start, flush=True)
            if record['bytes_downloaded'] != EXPECTED_BYTES:
                raise ValueError('Truncated movie')
            with partial.open('rb') as f:
                if f.read(1024) != header.read_bytes():
                    raise ValueError('Movie header changed')
            partial.rename(path)
            record.update(complete=True, sha256=h.hexdigest(), path=str(path.relative_to(ROOT)), seconds=time.perf_counter()-start)
            save(); print('complete', record['seconds'], record['sha256'], flush=True)
    except Exception as exc:
        record.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
