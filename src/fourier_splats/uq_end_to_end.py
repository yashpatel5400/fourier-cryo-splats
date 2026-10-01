"""Bookkeeping for paired local-alignment interval experiments."""
import numpy as np
from scipy.stats import beta


def dephase_observations(observations, q, shifts_A, field_A):
    y, q, shifts = (np.asarray(v, float) for v in (observations, q, shifts_A))
    if q.ndim != 3 or q.shape[-1] != 2 or y.shape != (len(q), 2*q.shape[1]) or shifts.shape != (len(q), 2):
        raise ValueError('Matching real/imaginary observations and shifts required')
    if not np.isfinite(field_A) or field_A <= 0 or not all(np.isfinite(v).all() for v in (y, q, shifts)):
        raise ValueError('Finite observations and positive field size required')
    nq = q.shape[1]
    value = (y[:, :nq]+1j*y[:, nq:])*np.exp(2j*np.pi*np.einsum('nqa,na->nq', q, shifts)/field_A)
    return np.concatenate([value.real, value.imag], axis=1).ravel()


def observed_interval(center, half_width, pilot_center, no_data_half_width, truth, fallback):
    values = [center, half_width, pilot_center, no_data_half_width, truth]
    if not np.isfinite(values).all() or min(half_width, no_data_half_width) < 0:
        raise ValueError('Finite interval parameters required')
    used = bool(fallback and half_width > no_data_half_width)
    selected_center = pilot_center if used else center
    selected_half = no_data_half_width if used else half_width
    return dict(raw_center=float(center), raw_half_width=float(half_width),
        raw_covered=bool(abs(center-truth) <= half_width),
        raw_correct_sign_exclusion=bool(truth != 0 and np.sign(truth)*center > half_width),
        center=float(selected_center), half_width=float(selected_half), fallback=used,
        covered=bool(abs(selected_center-truth) <= selected_half),
        correct_sign_exclusion=bool(truth != 0 and np.sign(truth)*selected_center > selected_half),
        relative_half_width=float(selected_half/no_data_half_width) if no_data_half_width else None,
        raw_relative_half_width=float(half_width/no_data_half_width) if no_data_half_width else None,
        half_width_over_abs_center=float(selected_half/abs(selected_center)) if selected_center else None,
        true_target=float(truth))


def binomial_interval(successes, trials, alpha=.05):
    if not isinstance(successes, (int, np.integer)) or not isinstance(trials, (int, np.integer)) or trials < 1 or not 0 <= successes <= trials or not 0 < alpha < 1:
        raise ValueError('Integer binomial counts and valid error budget required')
    return [0. if successes == 0 else float(beta.ppf(alpha/2, successes, trials-successes+1)),
            1. if successes == trials else float(beta.ppf(1-alpha/2, successes+1, trials-successes))]
