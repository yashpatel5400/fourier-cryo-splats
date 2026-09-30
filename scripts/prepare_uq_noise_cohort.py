#!/usr/bin/env python3
"""Reserve independent calibration exposures using source metadata only."""
import csv
import hashlib
import json
import pickle
from pathlib import Path
import numpy as np
from prepare_uq_splits import groups_for

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'research/uncertainty/confirmation/noise-calibration-v1'
COUNT = 128
SEED = 630929


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def select_representatives(groups, eligible, count, seed):
    """Uniform groups without replacement, then one uniform eligible row each."""
    groups = np.asarray(groups); eligible = np.asarray(eligible, bool)
    if groups.ndim != 1 or eligible.shape != groups.shape or count < 1:
        raise ValueError('Compatible group/eligibility vectors and positive count required')
    available = np.unique(groups[eligible])
    if len(available) < count:
        raise ValueError('Insufficient unused exposures; do not replace with seen groups')
    rng = np.random.default_rng(seed)
    chosen = rng.permutation(available)[:count]
    indices = np.array([rng.choice(np.flatnonzero(eligible & (groups == group)))
                        for group in chosen], dtype=np.int64)
    return np.sort(indices)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if any(OUT.glob('*-selection.csv')):
        raise RuntimeError('Preserve the existing cohort selection')
    for dataset in ['10028', '10049', '10076']:
        inputs = ROOT/f'background/cryodrgn_empiar/empiar{dataset}/inputs'
        cs_path = next(inputs.glob('*.cs')); cs = np.load(cs_path, allow_pickle=False)
        _, groups, detail = groups_for(dataset, np.arange(len(cs)))
        paths = [ROOT/f'research/uncertainty/splits/{dataset}.csv',
                 ROOT/f'research/uncertainty/confirmation/prediction-v1/{dataset}-selection.csv']
        excluded = set()
        for path in paths:
            with path.open() as file:
                excluded.update(row['source_group'] for row in csv.DictReader(file))
        eligible = ~np.isin(groups, sorted(excluded))
        filter_path = inputs/'filtered.ind.pkl'
        if filter_path.exists():
            published = np.zeros(len(cs), bool)
            with filter_path.open('rb') as file:
                published[np.asarray(pickle.load(file), dtype=int)] = True
            eligible &= published
        exclusions_path = ROOT/f'provenance/{dataset}-source-exclusions.json'
        exclusions = json.loads(exclusions_path.read_text()) if exclusions_path.exists() else []
        for exclusion in exclusions:
            eligible &= ~np.array([p.decode().endswith(exclusion['metadata_path_suffix']) for p in cs['blob/path']])
        selected = select_representatives(groups, eligible, COUNT, SEED+int(dataset))
        if len(np.unique(groups[selected])) != COUNT or set(groups[selected]) & excluded:
            raise AssertionError('Calibration groups must be distinct and previously unused')
        np.save(OUT/f'{dataset}-indices.npy', selected)
        path = OUT/f'{dataset}-selection.csv'
        with path.open('w') as file:
            writer = csv.writer(file, lineterminator='\n')
            writer.writerow(['output_index', 'metadata_source_index', 'source_group'])
            writer.writerows(zip(range(COUNT), selected, groups[selected]))
        manifest = {'stage': 'metadata-only reservation; no new particle pixels read',
            'dataset': dataset, 'count': COUNT, 'seed': SEED+int(dataset),
            'selection_policy': 'uniform random eligible exposure groups without replacement, then one uniform eligible particle per group; sort source indices',
            'eligible_groups': len(np.unique(groups[eligible])), 'eligible_particles': int(eligible.sum()),
            'excluded_groups': len(excluded), 'selected_groups': COUNT, 'excluded_group_overlap': 0,
            'selection_csv_sha256': sha(path), 'indices_npy_sha256': sha(OUT/f'{dataset}-indices.npy'),
            'excluded_selection_hashes': {str(p.relative_to(ROOT)): sha(p) for p in paths},
            'metadata_sha256': sha(cs_path), 'published_filter_sha256': sha(filter_path) if filter_path.exists() else None,
            'source_groups': detail, 'source_exclusions': exclusions,
            'source_code_sha256': sha(Path(__file__)),
            'release_condition': 'Do not download until all twelve weights and evaluation code are frozen and committed.'}
        (OUT/f'{dataset}-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        print(dataset, 'remaining groups', manifest['eligible_groups'], 'reserved', COUNT, 'overlap', 0, flush=True)


if __name__ == '__main__':
    main()
