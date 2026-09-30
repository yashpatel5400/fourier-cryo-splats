"""Scientific cohort selection and model-lock drift checks."""
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def script(name):
    import sys
    if str(ROOT/'scripts') not in sys.path:
        sys.path.insert(0, str(ROOT/'scripts'))
    spec = importlib.util.spec_from_file_location(name, ROOT/'scripts'/f'{name}.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_cohort_uses_groups_not_particle_multiplicities():
    module = script('prepare_uq_noise_cohort')
    groups = np.array(['excluded']*8+['a']*100+['b']+['c']*3)
    eligible = groups != 'excluded'
    selected = module.select_representatives(groups, eligible, 3, 7)
    assert set(groups[selected]) == {'a', 'b', 'c'}
    assert np.all(np.diff(selected) > 0)
    np.testing.assert_array_equal(selected, module.select_representatives(groups, eligible, 3, 7))
    with pytest.raises(ValueError, match='Insufficient'):
        module.select_representatives(groups, eligible, 4, 7)


def test_reserved_cohort_hashes_and_exposure_disjointness():
    module = script('prepare_uq_noise_cohort')
    import csv
    for dataset in ['10028', '10049', '10076']:
        path = module.OUT/f'{dataset}-selection.csv'
        rows = list(csv.DictReader(path.open()))
        manifest = json.loads((module.OUT/f'{dataset}-manifest.json').read_text())
        assert len(rows) == 128
        assert len({r['source_group'] for r in rows}) == 128
        assert module.sha(path) == manifest['selection_csv_sha256']
        excluded = set()
        for name, expected in manifest['excluded_selection_hashes'].items():
            prior = ROOT/name
            assert module.sha(prior) == expected
            excluded.update(r['source_group'] for r in csv.DictReader(prior.open()))
        assert not {r['source_group'] for r in rows} & excluded
        # The CSV is authoritative and reconstructs the ignored/numeric array.
        np.testing.assert_array_equal(np.load(module.OUT/f'{dataset}-indices.npy'),
                                      [int(r['metadata_source_index']) for r in rows])


def test_fresh_calibration_lock_detects_changed_input(tmp_path, monkeypatch):
    module = script('freeze_uq_noise_models')
    monkeypatch.setattr(module, 'ROOT', tmp_path); monkeypatch.setattr(module, 'COHORT', tmp_path)
    input_file = tmp_path/'weights'; input_file.write_bytes(b'locked scientific input')
    lock = {'target_lock_sha256': module.LOCK_SHA256, 'estimators': [{}]*12,
            'files': {'weights': hashlib.sha256(input_file.read_bytes()).hexdigest()}}
    (tmp_path/'locked-models.json').write_text(json.dumps(lock))
    module.verify_lock(require_committed=False)
    input_file.write_bytes(b'adapted using calibration data')
    with pytest.raises(ValueError, match='Frozen input changed'):
        module.verify_lock(require_committed=False)
