"""Explicit-seed entry point for future designs; the frozen study API is unchanged."""
from numbers import Integral
from .uq_cubic_optimization import optimize_cubic_weights


def fit_with_fresh_certificate(objective, initial_weights, *, certificate_seed, **options):
    """Require an explicit seed; callers must also prevent reuse across audits.

    A distinct integer is reproducibility hygiene, not a proof of statistical
    independence. The chosen matrix must not depend on its certificate probes.
    The frozen historical implementation retains its old signature for replay.
    """
    if isinstance(certificate_seed, bool) or not isinstance(certificate_seed, Integral) or certificate_seed < 0:
        raise ValueError('An explicit nonnegative integer certificate seed is required')
    if certificate_seed == options.get('optimization_seed', 650011):
        raise ValueError('Separate design and final-certificate random streams required')
    return optimize_cubic_weights(objective, initial_weights, certificate_seed=int(certificate_seed), **options)
