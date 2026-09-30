#!/usr/bin/env python3
"""Wait for the declared models, publish their lock, then run fresh calibration.

This orchestrates an already committed protocol. No feature, fit criterion,
cohort, statistical budget or inference code is selected using fresh outcomes.
"""
import argparse
import fcntl
import hashlib
import json
import subprocess
import sys
import time

from fourier_splats.uq_pilot_targets import read_target_lock
from freeze_uq_noise_models import ROOT, COHORT, BASE, sha, verify_lock


def invoke(*command):
    print('COMMAND', json.dumps(list(map(str, command))), flush=True)
    subprocess.run(list(map(str, command)), cwd=ROOT, check=True)


def complete_family(paths, deadline):
    last = None
    while True:
        completed = 0
        for path in paths:
            try:
                row = json.loads(path.read_text())
            except (FileNotFoundError, json.JSONDecodeError):
                continue
            if row.get('error'):
                raise RuntimeError(f'Failed source retained; inspect before continuing: {path}')
            completed += bool(row.get('complete'))
        if completed != last:
            print('COMPLETED_SOURCES', completed, 'of', len(paths), str(paths[0].parent), flush=True)
            last = completed
        if completed == len(paths):
            return
        if time.monotonic() >= deadline:
            raise TimeoutError('Declared source family did not finish within the waiting budget')
        time.sleep(20)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wait-hours', type=float, default=12)
    args = parser.parse_args()
    if not 0 < args.wait_hours <= 48:
        raise ValueError('Waiting budget must be positive and at most 48 hours')
    guard_path = ROOT/'tmp/fresh-noise-coordinator.lock'
    guard_path.parent.mkdir(exist_ok=True)
    with guard_path.open('a') as guard:
        fcntl.flock(guard, fcntl.LOCK_EX | fcntl.LOCK_NB)
        specification = read_target_lock(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
        pairs = [(d['dataset'], f['name']) for d in specification['datasets'] for f in d['features']]
        fits = [BASE/f'pilot-selected-fixed-v1/{ds}-{name}.json' for ds, name in pairs]
        applications = [BASE/f'pilot-selected-experimental-v1/{ds}-{name}.json' for ds, name in pairs]
        deadline = time.monotonic()+args.wait_hours*3600
        complete_family(fits, deadline)
        lock_path = COHORT/'locked-models.json'
        if not lock_path.exists():
            # Do not incorporate somebody else's staged work into the lock commit.
            subprocess.run(['git', 'diff', '--cached', '--quiet'], cwd=ROOT, check=True)
            invoke(sys.executable, ROOT/'scripts/freeze_uq_noise_models.py')
            lock = verify_lock(require_committed=False)
            for name, expected in lock['files'].items():
                if name.startswith(('src/', 'scripts/', 'research/')):
                    committed = subprocess.check_output(['git', 'show', f'HEAD:{name}'], cwd=ROOT)
                    if hashlib.sha256(committed).hexdigest() != expected:
                        raise RuntimeError(f'Commit the exact protocol/source before fresh pixels: {name}')
            names = [str(p.relative_to(ROOT)) for p in [lock_path, *fits]]
            invoke('git', 'add', '--', *names)
            invoke('git', 'commit', '-m', 'Freeze all twelve estimators before fresh noise calibration pixels')
        verify_lock(require_committed=True)
        invoke('git', 'push')
        print('PUBLISHED_MODEL_LOCK', sha(lock_path), flush=True)
        invoke(sys.executable, ROOT/'scripts/download_uq_noise_cohort.py')
        complete_family(applications, deadline)
        summary = ROOT/'results/uncertainty/confirmation/noise-calibration-v1/summary.json'
        if summary.exists():
            raise RuntimeError('Existing final or partial outcomes require inspection; never overwrite')
        invoke(sys.executable, ROOT/'scripts/apply_uq_fresh_noise.py')
        print('FRESH_CALIBRATION_CHECK_COMPLETE', flush=True)


if __name__ == '__main__':
    main()
