#!/usr/bin/env python3
"""Run the complete locked target grid, preserving completed and failed records."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import time
from fourier_splats.uq_pilot_targets import read_target_lock

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def invoke(script, *arguments, log=None):
    command = [sys.executable, str(ROOT/'scripts'/script), *map(str, arguments)]
    print('COMMAND', json.dumps(command), flush=True)
    subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT if log else None, check=True)


def available(path):
    if not path.exists():
        return False
    d = json.loads(path.read_text())
    if not d.get('complete') or d.get('error'):
        raise RuntimeError(f'Existing incomplete/failed case must be preserved and inspected: {path}')
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wait-first-fit', action='store_true')
    args = parser.parse_args()
    first = BASE/'pilot-selected-fixed-v1/10049-pilot_region_1.json'
    if args.wait_first_fit:
        print('WAIT_FOR_FIRST_COMPUTE_PILOT', flush=True)
        deadline = time.monotonic()+7200
        while time.monotonic() < deadline:
            try:
                d = json.loads(first.read_text()) if first.exists() else {}
            except json.JSONDecodeError:
                d = {}
            if d.get('error'):
                raise RuntimeError(f'First compute pilot failed: {d["error"]}')
            if d.get('complete'):
                print('FIRST_COMPUTE_PILOT_SECONDS', d['seconds'], flush=True)
                break
            time.sleep(10)
        else:
            raise TimeoutError('First compute pilot not complete after two hours; no grid launched')
    if not available(first):
        raise RuntimeError('Time the prescribed first case before launching the full grid')
    lock = read_target_lock(ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json')
    def dataset_study(di, dataset):
        ds = dataset['dataset']; log_path = ROOT/f'logs/uncertainty/pilot-selected-study-{ds}.log'
        with log_path.open('a') as log:
            invoke('run_uq_pilot_targets.py', '--datasets', ds, '--resume', log=log)
            for fi, feature in enumerate(dataset['features']):
                name = feature['name']
                for ai, angle in enumerate([0, 1, 2]):
                    path = BASE/'pilot-selected-pose-v1'/f'{ds}-{name}-{angle}.json'
                    if not available(path):
                        invoke('audit_uq_high_band_pose.py', '--fit', f'results/uncertainty/development/pilot-selected-fixed-v1/{ds}-{name}.json',
                               '--target', name, '--angle', angle, '--shift-A', .5, '--order', 80, '--design-order', 12,
                               '--iterations', 40, '--threads', 2, '--pilot-pairing', '--alpha-total', .05/12,
                               '--delta', 1e-6/12, '--certificate-seed', 620000+100*di+10*fi+ai,
                               '--output', 'pilot-selected-pose-v1', log=log)
                    refined = BASE/'pilot-selected-enclosing-v1'/path.name
                    if not available(refined):
                        invoke('probe_uq_ball_remainder.py', '--audit', str(path.relative_to(ROOT)),
                               '--domain', 'cube', '--output', 'pilot-selected-enclosing-v1', log=log)
                applied = BASE/'pilot-selected-experimental-v1'/f'{ds}-{name}.json'
                if not available(applied):
                    invoke('apply_uq_pilot_targets.py', '--datasets', ds, '--features', name, log=log)
        return ds
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(dataset_study, di, dataset) for di, dataset in enumerate(lock['datasets'])]
        errors = []
        for future in futures:
            try:
                print('COMPLETE_DATASET', future.result(), flush=True)
            except Exception as exc:
                errors.append(repr(exc))
                print('DATASET_FAILURE_RETAINED', repr(exc), flush=True)
        if errors:
            raise RuntimeError(f'Incomplete study; all independent datasets attempted: {errors}')
    invoke('summarize_uq_pilot_targets.py')
    print('COMPLETE_LOCKED_TARGET_GRID', flush=True)


if __name__ == '__main__':
    main()
