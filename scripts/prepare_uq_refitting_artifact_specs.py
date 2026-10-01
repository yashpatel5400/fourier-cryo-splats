#!/usr/bin/env python3
"""Hash-lock every prescribed local-refinement record, array, and frozen source."""
import argparse
import hashlib
import json
from pathlib import Path
from review_uq_candidate import snapshot_references
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--datasets',default='10028,10049,10076');args=ap.parse_args()
    for ds in args.datasets.split(','):
        if ds not in ['10028','10049','10076']:raise ValueError('Unknown geometry')
        owners={};files={}
        for study,count in [('local-alignment-calibration-v1',128),('end-to-end-local-pose-v1',200)]:
            directory=BASE/study/ds;owner=directory/'summary.json';d=json.loads(owner.read_text())
            if not d.get('complete') or len(d['records'])!=count:
                raise ValueError('Every prescribed trial must be attempted')
            owner_name=str(owner.relative_to(ROOT));owner_hash=sha(owner);owners[owner_name]=owner_hash
            # All partial arrays, exceptions and nonconverged outcomes are included.
            # The owner denotes a completed batch, not scientific success of every trial.
            for path in sorted(directory.iterdir()):
                if not path.is_file() or path.suffix not in ['.json','.npz']:
                    raise ValueError('Unexpected result file; explicitly review artifact scope')
                files[str(path.relative_to(ROOT))]=(path,owner_name,owner_hash)
                if path.suffix=='.json':
                    for ref in snapshot_references(json.loads(path.read_text())):
                        p=ROOT/ref
                        if sha(p)!=p.stem:raise ValueError('Archived source hash changed')
                        files[ref]=(p,owner_name,owner_hash)
            for entry in d['records']:
                rp=directory/f"replicate-{entry['replicate']:03d}.json"
                if sha(rp)!=entry['record_sha256']:raise ValueError('Changed replicate')
                rec=json.loads(rp.read_text());suffix='' if rec['complete'] else '-partial'
                arr=directory/f"replicate-{entry['replicate']:03d}{suffix}.npz"
                key='arrays_sha256' if rec['complete'] else 'partial_arrays_sha256'
                if sha(arr)!=rec[key]:raise ValueError('Changed replicate arrays')
        summary=BASE/'end-to-end-local-pose-summary-v1'/f'{ds}.json';d=json.loads(summary.read_text())
        if not d.get('complete') or d['attempted_replicates']!=200:raise ValueError('Verified final summary required')
        owner_name=str(summary.relative_to(ROOT));owner_hash=sha(summary)
        for path in [summary,summary.with_suffix('.csv')]:files[str(path.relative_to(ROOT))]=(path,owner_name,owner_hash)
        spec=dict(version=f'v0.7.0-dev-local-refitting-{ds}',files=[dict(path=name,sha256=sha(p),owner=owner,owner_sha256=oh)
            for name,(p,owner,oh) in sorted(files.items())],required_prior_arrays=[],
            scope='Complete 128-calibration and 200-test fixed-generator datasets for one acquisition geometry. All raw observations, poses, weights, intervals, failures and source snapshots retained. Completed batch owners do not imply every solver succeeded. This is local alignment, not global density learning or experimental truth coverage.',
            reproduction='See END-TO-END-LOCAL-POSE-PROTOCOL.md, frozen source snapshots and per-replicate seeds. Reproduce scientific fits from generator arrays; summary script verifies every array/record hash.')
        path=ROOT/f'provenance/uncertainty/refitting-{ds}-artifact-spec-v1.json'
        if path.exists():raise RuntimeError('Preserve previous artifact specification')
        path.write_text(json.dumps(spec,indent=2)+'\n');print(ds,len(files),'artifact members',flush=True)


if __name__=='__main__':main()
