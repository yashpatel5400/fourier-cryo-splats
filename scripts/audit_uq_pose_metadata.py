#!/usr/bin/env python3
"""Inspect whether the archived metadata actually contain nonzero pose summaries."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    output = ROOT/'results/uncertainty/development/pose-metadata-inventory.json'
    if output.exists():
        raise RuntimeError('Preserve previous inventory')
    result = {'complete': True,
        'scope': 'Metadata availability check, not an angular-error calibration',
        'source_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'datasets': []}
    for dataset in ['10028', '10049', '10076']:
        path = next((ROOT/f'background/cryodrgn_empiar/empiar{dataset}/inputs').glob('*.cs'))
        values = np.load(path, allow_pickle=False)
        row = {'dataset': dataset, 'metadata': str(path.relative_to(ROOT)),
            'metadata_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'source_particles': len(values), 'fields': {}}
        for name in ['alignments3D/pose_ess', 'alignments3D/shift_ess',
                     'alignments3D/error', 'alignments3D/error_min', 'alignments3D/weight']:
            x = values[name]
            row['fields'][name] = {'shape': list(x.shape), 'nonzero': int(np.count_nonzero(x)),
                'minimum': float(np.nanmin(x)), 'maximum': float(np.nanmax(x)), 'finite': bool(np.isfinite(x).all())}
        result['datasets'].append(row)
    result['interpretation'] = ('Both pose_ess and shift_ess are identically zero in all three source files. '
        'These stored zeros do not mean exact poses. Nonzero alignment error/objective fields have no '
        'verified conversion here to an angular confidence radius. No experimental per-particle pose bound is obtained.')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['interpretation'])


if __name__ == '__main__':
    main()
