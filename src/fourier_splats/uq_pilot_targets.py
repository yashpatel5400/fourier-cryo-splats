"""Validation of the pilot-only target lock, with no reference-image access."""
import hashlib
import json
import numpy as np
from .uq_continuous import cell_target_coefficients


LOCK_SHA256 = '57122bb2ee0d545584d854c1e8845adeced7620d4926a923fa9bf2c74af3e834'


def read_target_lock(path):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != LOCK_SHA256:
        raise ValueError('Pilot target lock changed')
    lock = json.loads(raw)
    if [d['dataset'] for d in lock['datasets']] != ['10028', '10049', '10076']:
        raise ValueError('Unexpected target family')
    for row in lock['datasets']:
        validate_target_row(row)
    return lock


def validate_target_row(row):
    names = ['pilot_region_1', 'pilot_region_2', 'pilot_region_3', 'matched_center']
    if row['selected_count'] != 3 or [f['name'] for f in row['features']] != names:
        raise ValueError('All three pilot targets and the center control are required')
    field = float(row['field_A'])
    if not np.isfinite(field) or field <= 0 or row['width_A'] != 20.:
        raise ValueError('Invalid physical scale')
    np.testing.assert_allclose(row['width_fraction_field']*field, 20., rtol=1e-13)
    centers = np.array([f['center_fraction_field'] for f in row['features']])
    if centers.shape != (4, 3) or not np.isfinite(centers).all():
        raise ValueError('Invalid target coordinates')
    if np.max(np.linalg.norm(centers, axis=1)) > .30 or np.any(centers[-1] != 0):
        raise ValueError('Target outside selection region or misplaced center control')
    for i in range(3):
        for j in range(i):
            if np.linalg.norm(centers[i]-centers[j])*field < 40.-1e-10:
                raise ValueError('Pilot target separation is smaller than specified')
    for f in row['features']:
        np.testing.assert_allclose(np.asarray(f['center_fraction_field'])*field, f['center_A'], atol=1e-12)


def verify_pilot_scores(row, pilot):
    if np.shape(pilot) != (24**3,):
        raise ValueError('The locked pilot is the full 24-cell function')
    np.testing.assert_allclose(np.linalg.norm(pilot), 1., rtol=1e-12)
    for feature in row['features']:
        value = cell_target_coefficients(24, [feature['center_fraction_field']], [1],
                                         row['width_fraction_field'])@pilot
        np.testing.assert_allclose(value, feature['pilot_expected_feature'], rtol=2e-12, atol=1e-12)
