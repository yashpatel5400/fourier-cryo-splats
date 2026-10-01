#!/usr/bin/env python3
"""Audit the existing acquisition joins; preserve original analysis and splits."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from prepare_uq_splits import groups_for

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/uncertainty/development/background-source-group-correction-v1'


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    if OUT.exists():
        raise ValueError('Preserve prior correction audits')
    OUT.mkdir(parents=True)
    old = json.loads((ROOT/'research/uncertainty/splits/manifest.json').read_text())
    result = {'complete': False, 'datasets': [], 'purpose': 'Metadata correction; no spectrum recomputation or statistical claim'}
    for ds in ['10028', '10049', '10076']:
        selected, groups, detail = groups_for(ds)
        path = ROOT / 'research/uncertainty/splits' / (ds + '.csv')
        with path.open() as f:
            rows = list(csv.DictReader(f))
        assert sha(path) == old['datasets'][ds]['split_file_sha256']
        np.testing.assert_array_equal([int(r['output_index']) for r in rows], np.arange(8192))
        np.testing.assert_array_equal([int(r['metadata_source_index']) for r in rows], selected)
        np.testing.assert_array_equal([r['source_group'] for r in rows], groups)
        unique, counts = np.unique(groups, return_counts=True)
        assert len(unique) == old['datasets'][ds]['selected_groups']
        # Check whole-group assignments from saved labels, independent of hashing routine.
        assignment = {}
        for g, r in zip(groups, rows):
            assignment.setdefault(g, set()).add(r['split'])
        assert all(len(v) == 1 for v in assignment.values())
        result['datasets'].append(dict(dataset=ds, particles=len(groups),
            selected_source_groups=len(unique), minimum_particles_per_group=int(counts.min()),
            median_particles_per_group=float(np.median(counts)), maximum_particles_per_group=int(counts.max()),
            source=detail, split_file_sha256=sha(path), all_joins_match_saved_splits=True,
            no_source_group_crosses_saved_split=True))
    result['complete'] = True
    (OUT/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
