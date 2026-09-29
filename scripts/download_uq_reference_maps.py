#!/usr/bin/env python3
"""Download deposited maps as known generators for semisynthetic experiments.

These maps are not ground truth for the original experimental particles. Their
role is to define a reproducible known signal for fresh synthetic observations.
"""
import concurrent.futures
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import mrcfile
import requests

ROOT=Path(__file__).resolve().parents[1]
REFERENCES={'10028':'2660','10049':'6487','10076':'8434'}


def download(item):
    empiar,emd=item
    destination=ROOT/'data/uncertainty/references'/('emd_'+emd+'.map')
    destination.parent.mkdir(parents=True,exist_ok=True)
    source=f'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{emd}/map/emd_{emd}.map.gz'
    compressed=destination.with_suffix('.map.gz')
    if not compressed.exists():
        with requests.get(source,stream=True,timeout=90) as response:
            response.raise_for_status()
            partial=compressed.with_suffix('.part')
            with partial.open('wb') as out:
                for block in response.iter_content(2**20):out.write(block)
            partial.replace(compressed)
    if not destination.exists():
        with gzip.open(compressed,'rb') as source_file,destination.open('wb') as out:
            while block:=source_file.read(2**20):out.write(block)
    def digest(path):
        h=hashlib.sha256()
        with path.open('rb') as f:
            while block:=f.read(2**20):h.update(block)
        return h.hexdigest()
    with mrcfile.open(destination,permissive=False) as m:
        record={'empiar':empiar,'emdb':emd,'source':source,
                'role':'known semisynthetic generator, not experimental ground truth',
                'shape_zyx':list(m.data.shape),'voxel_size_xyz_A':[float(m.voxel_size[x]) for x in ['x','y','z']],
                'axis_mapping':[int(m.header[x]) for x in ['mapc','mapr','maps']],
                'origin_xyz':[float(m.header.origin[x]) for x in ['x','y','z']],
                'start_xyz':[int(m.header[x]) for x in ['nxstart','nystart','nzstart']],
                'sha256_map':digest(destination),'sha256_gzip':digest(compressed),
                'downloaded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'bytes_map':destination.stat().st_size,'bytes_gzip':compressed.stat().st_size}
    print(emd,record['shape_zyx'],record['voxel_size_xyz_A'],flush=True)
    return record


if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(3) as pool:
        records=list(pool.map(download,REFERENCES.items()))
    target=ROOT/'research/uncertainty/reference-map-manifest.json'
    target.write_text(json.dumps(records,indent=2)+'\n')
