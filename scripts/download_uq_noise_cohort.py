#!/usr/bin/env python3
"""Download only after the complete estimator/cohort lock has been committed."""
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from freeze_uq_noise_models import ROOT, COHORT, sha, verify_lock


def main():
    verify_lock()
    lock_sha = sha(COHORT/'locked-models.json')
    outputs = ROOT/'provenance/uncertainty/noise-calibration-v1'
    outputs.mkdir(parents=True, exist_ok=True)
    for dataset in ['10028', '10049', '10076']:
        out = ROOT/f'data/uncertainty/confirmation/noise-calibration-v1/{dataset}'
        record = outputs/f'{dataset}-download.json'
        if record.exists():
            prior = json.loads(record.read_text())
            if prior['model_lock_sha256'] != lock_sha or any(sha(ROOT/n) != h for n, h in prior['cached_files'].items()):
                raise ValueError('Previously completed download changed')
            print('REUSE', dataset, flush=True); continue
        subprocess.run([sys.executable, str(ROOT/'scripts/download_data.py'), dataset,
            '--count', '128', '--box', '64', '--workers', '4', '--resume',
            '--indices-file', str(COHORT/f'{dataset}-indices.npy'), '--output-directory', str(out)], cwd=ROOT, check=True)
        np.testing.assert_array_equal(np.load(out/'indices.npy'), np.load(COHORT/f'{dataset}-indices.npy'))
        images = np.load(out/'images.npy', mmap_mode='r')
        if images.shape != (128, 64, 64) or not np.isfinite(images).all():
            raise ValueError('Incomplete or invalid cropped images')
        manifest = json.loads((out/'manifest.json').read_text())
        if manifest['selection_file_sha256'] != sha(COHORT/f'{dataset}-indices.npy'):
            raise ValueError('Downloaded the wrong cohort')
        provenance = {'complete': True, 'dataset': dataset, 'model_lock_sha256': lock_sha,
            'selection_csv_sha256': sha(COHORT/f'{dataset}-selection.csv'),
            'cached_files': {str((out/name).relative_to(ROOT)): sha(out/name)
                             for name in ['images.npy', 'indices.npy', 'metadata.npz', 'manifest.json']},
            'download_manifest': manifest}
        record.write_text(json.dumps(provenance, indent=2)+'\n')
        print('VERIFIED', dataset, '128 previously unused exposures', flush=True)


if __name__ == '__main__':
    main()
