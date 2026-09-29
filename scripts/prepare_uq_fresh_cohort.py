#!/usr/bin/env python3
"""Select additional exposure groups using only metadata, before image access."""
import csv
import hashlib
import json
import pickle
from pathlib import Path
import numpy as np
from prepare_uq_splits import groups_for
ROOT = Path(__file__).resolve().parents[1]
SEED = 609421
COUNT = 4096


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    out = ROOT/'research/uncertainty/confirmation/prediction-v1'
    out.mkdir(parents=True, exist_ok=True)
    for dataset in ['10028', '10049', '10076']:
        inputs = ROOT/'background/cryodrgn_empiar'/('empiar'+dataset)/'inputs'
        cs_path = next(inputs.glob('*.cs')); cs = np.load(cs_path, allow_pickle=False)
        _, groups, detail = groups_for(dataset, np.arange(len(cs)))
        previous_path = ROOT/'research/uncertainty/splits'/f'{dataset}.csv'
        previous = list(csv.DictReader(previous_path.open()))
        excluded_groups = set(row['source_group'] for row in previous)
        eligible = ~np.isin(groups, sorted(excluded_groups))
        if (inputs/'filtered.ind.pkl').exists():
            filtered = np.zeros(len(cs), dtype=bool)
            filtered[np.asarray(pickle.load((inputs/'filtered.ind.pkl').open('rb')), dtype=int)] = True
            eligible &= filtered
        exclusions_path = ROOT/'provenance'/f'{dataset}-source-exclusions.json'
        exclusions = json.loads(exclusions_path.read_text()) if exclusions_path.exists() else []
        for exclusion in exclusions:
            eligible &= ~np.array([p.decode().endswith(exclusion['metadata_path_suffix']) for p in cs['blob/path']])
        starts = np.arange(0, len(cs), 64); np.random.default_rng(SEED+int(dataset)).shuffle(starts)
        selected = []
        for start in starts:
            selected.extend(i for i in range(start, min(start+64, len(cs))) if eligible[i])
            if len(selected) >= COUNT:
                break
        selected = np.sort(np.array(selected[:COUNT], dtype=np.int64))
        if len(selected) != COUNT:
            raise ValueError('Insufficient new-source-group particles')
        if set(groups[selected]) & excluded_groups:
            raise AssertionError('Exposure-group leakage')
        # Numeric NPY is an execution convenience. The public CSV is authoritative
        # and can recreate the same array without redistributing any image bytes.
        np.save(out/f'{dataset}-indices.npy', selected)
        path = out/f'{dataset}-selection.csv'
        with path.open('w') as file:
            writer = csv.writer(file, lineterminator='\n')
            writer.writerow(['output_index', 'metadata_source_index', 'source_group'])
            writer.writerows(zip(range(len(selected)), selected, groups[selected]))
        manifest = {'stage': 'metadata-only selection before fresh image access or outcome evaluation',
                    'dataset': dataset, 'seed': SEED+int(dataset), 'count': COUNT,
                    'selected_groups': len(np.unique(groups[selected])), 'excluded_development_groups': len(excluded_groups),
                    'eligible_new_group_particles': int(eligible.sum()), 'group_overlap': 0,
                    'policy': 'random permutation of source blocks of 64; published filter and source exclusions; omit every previously seen source group; retain first 4096 eligible particles, sort source indices',
                    'selection_csv_sha256': digest(path), 'indices_npy_sha256': digest(out/f'{dataset}-indices.npy'),
                    'development_split_sha256': digest(previous_path), 'source_metadata_sha256': digest(cs_path),
                    'source_groups': detail, 'source_exclusions': exclusions}
        (out/f'{dataset}-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        print(dataset, 'particles', COUNT, 'fresh groups', manifest['selected_groups'], 'overlap', manifest['group_overlap'], flush=True)


if __name__ == '__main__':
    main()
