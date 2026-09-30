import copy
from pathlib import Path
import numpy as np
import pytest
from fourier_splats.uq_pilot_targets import read_target_lock, validate_target_row, verify_pilot_scores


ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json'


def test_lock_and_physical_separation():
    lock = read_target_lock(LOCK)
    assert sum(len(d['features']) for d in lock['datasets']) == 12
    for row in lock['datasets']:
        xyz = np.array([f['center_A'] for f in row['features'][:3]])
        assert min(np.linalg.norm(xyz[i]-xyz[j]) for i in range(3) for j in range(i)) >= 40.-1e-10


def test_lock_rejects_changes(tmp_path):
    path = tmp_path/'lock.json'; path.write_bytes(LOCK.read_bytes()+b' ')
    with pytest.raises(ValueError, match='changed'):
        read_target_lock(path)


def test_reject_colliding_targets():
    row = copy.deepcopy(read_target_lock(LOCK)['datasets'][0])
    row['features'][1]['center_fraction_field'] = row['features'][0]['center_fraction_field']
    with pytest.raises(ValueError, match='separation'):
        validate_target_row(row)


def test_independent_score_replay():
    # Replay all locked scores from the original pilot inputs, with no map read.
    import sys
    sys.path.insert(0, str(ROOT/'scripts'))
    from audit_uq_grid_refinement import model
    from fourier_splats.uq_data import particle_geometry
    for row in read_target_lock(LOCK)['datasets']:
        dataset = row['dataset']
        g = particle_geometry(ROOT, dataset, 'inference_half0', radius=5, count=128, seed=609315)
        cp = ROOT/'results/uncertainty/development/representation'/dataset/'real_particles-spacing-2.0.npz'
        if not cp.exists():
            pytest.skip('Archived pilot arrays are needed for integration replay')
        op, pilot, _, _ = model(g, np.load(cp), 24)
        verify_pilot_scores(row, op.expand(pilot))
