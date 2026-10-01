#!/usr/bin/env python3
"""Independently recompute every saved interval decision from release records.

This arithmetic audit does not invoke the fitting or reporting modules and
makes no statement about whether the assumed class or simulated truths are
scientifically adequate. The signed endpoints are the scientific quantities
being checked, not an inference from a saved coverage Boolean.
"""
from pathlib import Path
import hashlib
import json
import math
import tarfile
import argparse


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise ValueError('Preserve prior verification')
    root=Path(__file__).resolve().parents[1];records=[]
    for ds in ['10028','10049','10076']:
        archive=root/f'output/artifacts/v0.7.0-dev-refitting-{ds}.tar.gz'
        expected={f'results/uncertainty/development/end-to-end-local-pose-v1/{ds}/replicate-{i:03d}.json' for i in range(200)}
        cases=0;seen=set();max_margin=-math.inf;min_margin=math.inf
        with tarfile.open(archive,'r|gz') as stream:
            for member in stream:
                if member.name not in expected:continue
                assert member.name not in seen;seen.add(member.name)
                with stream.extractfile(member) as f:d=json.load(f)
                assert d['complete'] and len(d['intervals'])==72
                cells=set()
                for r in d['intervals']:
                    key=tuple(r[x] for x in ['template','target','method','image_mode'])
                    assert key not in cells;cells.add(key)
                    truth=r['true_target'];assert math.isfinite(truth)
                    for prefix in ['', 'raw_']:
                        c=r[prefix+'center'];w=r[prefix+'half_width']
                        assert math.isfinite(c) and math.isfinite(w) and w>=0
                        covered=(c-w<=truth<=c+w)
                        assert covered==r[prefix+'covered'],(ds,d['replicate'],key,prefix)
                        signed=((truth>0 and c-w>0) or (truth<0 and c+w<0))
                        assert signed==r[prefix+'correct_sign_exclusion'],(ds,d['replicate'],key,prefix)
                        margin=w-abs(c-truth);min_margin=min(min_margin,margin);max_margin=max(max_margin,margin)
                    cases+=1
        assert seen==expected and cases==14400
        records.append(dict(dataset=ds,replicates=200,intervals=cases,raw_and_selected_decisions=2*cases,
            all_saved_coverage_and_sign_flags_match=True,min_absolute_coverage_margin=min_margin,
            max_absolute_coverage_margin=max_margin))
    result=dict(complete=True,scope='Independent arithmetic check of every archived raw/selected interval and correct-sign flag. No fitting, no new statistical validity claim, no alteration of the frozen results.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),records=records)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
