#!/usr/bin/env python3
"""Finish current local jobs: verify, summarize and package each complete batch.

This watches only the three already launched studies. It does not submit jobs,
change scientific parameters, retry failed trials, upload anything or declare
research success. Every wait is bounded and all subprocess failures stop it.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
DATASETS=['10028','10049','10076']


def run(*args):
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    subprocess.run([sys.executable,*args],cwd=ROOT,env=env,check=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--max-wait-seconds',type=int,default=7200);args=parser.parse_args()
    if not 0<=args.max_wait_seconds<=14400:raise ValueError('Bounded watch duration required')
    start=time.monotonic();finished=set()
    while len(finished)<3:
        for ds in DATASETS:
            if ds in finished:continue
            p=ROOT/f'results/uncertainty/development/end-to-end-local-pose-v1/{ds}/summary.json'
            try:d=json.loads(p.read_text())
            except json.JSONDecodeError:continue  # Writer may be replacing a progress record.
            if not d.get('complete'):continue
            if len(d['records'])!=200:raise ValueError('Completed flag without all 200 attempts')
            # These scripts exclusively create their final outputs. Existing output
            # is not silently overwritten or reused in a different logical run.
            run('scripts/summarize_uq_end_to_end_local.py','--datasets',ds)
            run('scripts/prepare_uq_refitting_artifact_specs.py','--datasets',ds)
            run('scripts/package_uq_increment.py',f'provenance/uncertainty/refitting-{ds}-artifact-spec-v1.json',
                '--output',f'output/artifacts/v0.7.0-dev-refitting-{ds}.tar.gz')
            finished.add(ds);print('FINALIZED',ds,flush=True)
        if len(finished)==3:break
        if time.monotonic()-start>=args.max_wait_seconds:raise TimeoutError('Studies still running; all existing results preserved')
        time.sleep(30)
    run('scripts/write_uq_end_to_end_figures.py')
    print('All three local-refinement batches verified, packaged and reported. Scientific acceptance remains a separate assessment.',flush=True)


if __name__=='__main__':main()
