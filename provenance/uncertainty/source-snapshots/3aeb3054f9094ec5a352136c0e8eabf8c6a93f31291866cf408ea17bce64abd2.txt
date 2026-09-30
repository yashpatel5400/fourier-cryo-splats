#!/usr/bin/env python3
"""Freeze the whole twelve-estimator family before fresh calibration access."""
import hashlib
import json
from pathlib import Path
import subprocess
from fourier_splats.uq_pilot_targets import read_target_lock, LOCK_SHA256

ROOT = Path(__file__).resolve().parents[1]
COHORT = ROOT/'research/uncertainty/confirmation/noise-calibration-v1'
BASE = ROOT/'results/uncertainty/development'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_lock(require_committed=True):
    path = COHORT/'locked-models.json'
    lock = json.loads(path.read_text())
    if lock['target_lock_sha256'] != LOCK_SHA256 or len(lock['estimators']) != 12:
        raise ValueError('Wrong or incomplete estimator family')
    for name, expected in lock['files'].items():
        if sha(ROOT/name) != expected:
            raise ValueError(f'Frozen input changed: {name}')
    if require_committed:
        relative = str(path.relative_to(ROOT))
        published = subprocess.check_output(['git', 'show', f'HEAD:{relative}'], cwd=ROOT)
        if published != path.read_bytes():
            raise ValueError('Commit the final lock before accessing new pixels')
    return lock


def main():
    path = COHORT/'locked-models.json'
    if path.exists():
        raise RuntimeError('Preserve the existing lock')
    fresh = ROOT/'data/uncertainty/confirmation/noise-calibration-v1'
    if fresh.exists() and any(fresh.rglob('images.npy')):
        raise RuntimeError('Cannot freeze prospectively after new image access')
    target = ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json'
    specification = read_target_lock(target)
    files = [target, *sorted(COHORT.glob('*.md')), *sorted(COHORT.glob('*-selection.csv')),
             *sorted(COHORT.glob('*-manifest.json')), *sorted(COHORT.glob('*-indices.npy')),
             *sorted((ROOT/'src/fourier_splats').glob('*.py'))]
    scripts = ['freeze_uq_noise_models.py', 'prepare_uq_noise_cohort.py',
               'download_uq_noise_cohort.py', 'apply_uq_fresh_noise.py', 'download_data.py',
               'apply_uq_pilot_targets.py', 'run_uq_pilot_targets.py',
               'audit_uq_high_band_pose.py', 'probe_uq_ball_remainder.py', 'audit_uq_grid_refinement.py']
    files += [ROOT/'scripts'/name for name in scripts]
    estimators = []
    for dataset in specification['datasets']:
        ds = dataset['dataset']
        manifest = json.loads((COHORT/f'{ds}-manifest.json').read_text())
        if manifest['selected_groups'] != 128 or manifest['excluded_group_overlap'] != 0:
            raise ValueError('Cohort selection does not match protocol')
        for filename in ['images.npy', 'indices.npy', 'metadata.npz', 'manifest.json']:
            files.append(ROOT/f'data/{ds}/{filename}')
        files += [ROOT/f'research/uncertainty/splits/{ds}.csv',
                  ROOT/f'research/uncertainty/confirmation/prediction-v1/{ds}-selection.csv',
                  BASE/f'experimental-noise-grouped/{ds}.json',
                  BASE/f'representation/{ds}/real_particles-spacing-2.0.npz']
        for feature in dataset['features']:
            name = feature['name']; fp = BASE/f'pilot-selected-fixed-v1/{ds}-{name}.json'
            fit = json.loads(fp.read_text())
            if not fit.get('complete') or fit.get('error') or fit['target_lock_sha256'] != LOCK_SHA256:
                raise ValueError(f'Incomplete/mismatched estimator: {fp}')
            row = next(r for r in fit['targets'] if r['target'] == name)
            wp = fp.with_name(f'{ds}-{name}-{row["width_fraction_field"]}-weights.npz')
            if sha(wp) != row['source_weights_sha256']:
                raise ValueError('Estimator weights changed')
            files += [fp, wp]
            estimators.append({'dataset': ds, 'feature': name,
                'fit': str(fp.relative_to(ROOT)), 'weights': str(wp.relative_to(ROOT)),
                'relative_objective_gap': row['fit']['history'][-1]['relative_gap'],
                'converged': row['fit']['converged']})
    contents = {'stage': 'all twelve estimators frozen before fresh calibration image access',
        'target_lock_sha256': LOCK_SHA256, 'estimators': estimators,
        'git_head_before_lock': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'files': {str(p.relative_to(ROOT)): sha(p) for p in dict.fromkeys(files)},
        'scope': 'Fresh calibration only; inference images and nuisance assumptions are not newly validated.'}
    path.write_text(json.dumps(contents, indent=2)+'\n')
    verify_lock(require_committed=False)
    print('LOCKED', len(estimators), 'estimators;', len(contents['files']), 'files; commit before download', flush=True)


if __name__ == '__main__':
    main()
