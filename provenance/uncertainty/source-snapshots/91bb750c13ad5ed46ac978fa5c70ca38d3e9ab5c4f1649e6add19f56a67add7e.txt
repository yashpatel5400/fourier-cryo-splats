"""Finite fixed-design Gaussian contrasts with selection accounted for."""
import numpy as np
from scipy.linalg import null_space
from .uq_group_noise import grouped_estimator_variance_upper


def fixed_noise_contrasts(ctf_design, pilot_design, relative_tolerance=1e-10):
    ctf = np.asarray(ctf_design, float); pilot = np.asarray(pilot_design, float)
    if ctf.ndim != 2 or pilot.ndim != 2 or len(ctf) != len(pilot) or len(ctf) < 18:
        raise ValueError('Compatible fixed designs with at least eighteen rows required')
    if not np.isfinite(ctf).all() or not np.isfinite(pilot).all() or not 0 < relative_tolerance < 1:
        raise ValueError('Finite fixed designs and positive rank threshold required')
    n = len(ctf); constant = np.ones((n, 1))/np.sqrt(n)
    family = {'raw': np.eye(n), 'mean': null_space(constant.T).T}; diagnostics = {}
    for name, design, counts in [('ctf', ctf, [3, 7]), ('pilot', pilot, [3, 7, 15])]:
        centered = design-design.mean(axis=0)
        u, singular, _ = np.linalg.svd(centered, full_matrices=False)
        available = int(np.sum(singular > relative_tolerance*singular[0])) if len(singular) and singular[0] > 0 else 0
        for count in counts:
            retained = min(count, available); key = f'{name}_{count+1}'
            removed = np.column_stack([constant, u[:, :retained]])
            q = null_space(removed.T, rcond=relative_tolerance).T
            if len(q) != n-retained-1: raise FloatingPointError('Unexpected contrast rank')
            family[key] = q
            diagnostics[key] = dict(requested_removed_rank=count+1, removed_rank=retained+1,
                available_centered_rank=available, singular_values=singular.tolist(),
                relative_rank_tolerance=relative_tolerance,
                removed_design_residual=float(np.linalg.norm(q@removed)),
                orthonormality_error=float(np.linalg.norm(q@q.T-np.eye(len(q)))))
    for name,q in family.items():
        np.testing.assert_allclose(q@q.T,np.eye(len(q)),rtol=1e-11,atol=1e-11)
        if name not in diagnostics:
            diagnostics[name]=dict(removed_rank=n-len(q),orthonormality_error=float(np.linalg.norm(q@q.T-np.eye(len(q)))))
    return family, diagnostics


def projected_grouped_variance_upper(observations, raw_weights, groups, contrasts, failure_probability=.005):
    y = np.asarray(observations, float)
    if y.ndim != 2 or not np.isfinite(y).all() or not contrasts:
        raise ValueError('Finite calibration rows and nonempty fixed contrast family required')
    if not 0 < failure_probability < 1: raise ValueError('Invalid error allowance')
    records = []
    for name,q in contrasts.items():
        q=np.asarray(q,float)
        if q.ndim != 2 or q.shape[1]!=len(y) or len(q)<1: raise ValueError('Invalid contrast shape')
        np.testing.assert_allclose(q@q.T,np.eye(len(q)),rtol=1e-10,atol=1e-10)
        bound=grouped_estimator_variance_upper(q@y,raw_weights,groups,failure_probability/len(contrasts))
        records.append(dict(name=name,independent_contrast_rows=len(q),bound=bound))
    chosen=min(records,key=lambda r:r['bound']['estimator_variance_upper'])
    return dict(noise_sd_upper=chosen['bound']['noise_sd_upper'],
        estimator_variance_upper=chosen['bound']['estimator_variance_upper'],selected=chosen['name'],
        family_size=len(contrasts),family_failure_probability=failure_probability,
        individual_failure_probability=failure_probability/len(contrasts),records=records,
        assumptions='All contrasts/weights fixed independently of calibration/inference noise; independent Gaussian calibration rows with one common covariance and arbitrary means; independent inference exposures with arbitrary joint Gaussian dependence inside each exposure.')
