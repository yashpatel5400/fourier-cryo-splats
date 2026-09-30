#!/usr/bin/env python3
"""Complete the remaining bounded-grid diagnostics as each saved fit finishes."""
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'
FOLDER = 'pose-aware-optimized-shift05-bounded'


def invoke(script, *arguments):
    cmd = [sys.executable, str(ROOT/'scripts'/script), *map(str, arguments)]
    print('COMMAND', json.dumps(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def main():
    for angle in ['0.5', '1', '2']:
        name = f'10028-contrast-{angle}.json'; source = BASE/FOLDER/name
        deadline = time.monotonic()+10800
        while time.monotonic() < deadline:
            try:
                d = json.loads(source.read_text()) if source.exists() else {}
            except json.JSONDecodeError:
                d = {}
            if d.get('error'):
                raise RuntimeError(f'Fit failed; retain and inspect: {source}')
            if d.get('complete'):
                break
            time.sleep(10)
        else:
            raise TimeoutError(f'No completion within three hours: {source}')
        rel = str(source.relative_to(ROOT))
        stages = [('pose-optimized-diagnostics', 'evaluate_uq_optimized_pose.py', []),
                  ('joint-bias-audit', 'audit_uq_joint_bias.py', []),
                  ('joint-bias-sharp-audit', 'audit_uq_joint_bias.py', ['--sharp-cubic']),
                  ('joint-bias-pilot-sharp-audit', 'audit_uq_joint_bias.py', ['--sharp-cubic', '--pilot-pairing'])]
        for directory, script, arguments in stages:
            path = BASE/directory/FOLDER/name
            if path.exists():
                previous = json.loads(path.read_text())
                if not previous.get('complete') or previous.get('error'):
                    raise RuntimeError(f'Preserve incomplete prior audit: {path}')
                continue
            invoke(script, '--fit', rel, *arguments)
    invoke('summarize_uq_joint_bias.py')
    print('COMPLETE_REMAINING_BOUNDED_GRID_DIAGNOSTICS', flush=True)


if __name__ == '__main__':
    main()
