"""Known-pose population diagnostics, not a latent-pose uncertainty method."""
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logit, roots_hermitenorm


def amplitude_tangent_information(full, region):
    """Real-amplitude Schur complements at the full and deleted means."""
    full, region = np.asarray(full), np.asarray(region)
    if full.shape != region.shape or full.ndim != 2:
        raise ValueError('Matching image by Fourier-coordinate arrays required')
    energy = np.sum(abs(region) ** 2, axis=1)
    residuals = []
    for mean in [full, full - region]:
        norm = np.sum(abs(mean) ** 2, axis=1)
        cross = np.real(np.sum(np.conj(mean) * region, axis=1))
        projection = np.divide(cross ** 2, norm, out=np.zeros_like(norm), where=norm > 0)
        residual = energy - projection
        if np.any(residual < -1e-12 * np.maximum(1, energy)):
            raise ArithmeticError('Negative real-amplitude Schur complement')
        residuals.append(np.maximum(residual, 0))
    return energy, *residuals


@lru_cache(maxsize=4)
def normal_quadrature(nodes):
    z, w = roots_hermitenorm(nodes)
    return z, w / np.sqrt(2 * np.pi)


def gaussian_mixture_information(energy, fraction, nodes=64):
    """Exact-in-model J_t(s), with deterministic 1D numerical quadrature."""
    s = np.asarray(energy, float)
    if s.ndim != 1 or np.any(s < 0) or not np.isfinite(s).all() or not 0 < fraction < 1:
        raise ValueError('Nonnegative finite energies and interior fraction required')
    z, w = normal_quadrature(nodes)
    out = np.empty_like(s)
    for begin in range(0, len(s), 4096):
        ss = s[begin:begin + 4096, None]
        base = np.sqrt(ss) * z + logit(fraction)
        difference = expit(base + ss / 2) - expit(base - ss / 2)
        out[begin:begin + len(ss)] = (difference @ w) / (fraction * (1 - fraction))
    return out


def pseudo_true_fraction(energy, law_a, law_b, fraction=.75):
    """Root of the expected pooled conditional score; independent 128-node check."""
    wa, wb = np.asarray(law_a, float), np.asarray(law_b, float)
    for w in [wa, wb]:
        if w.shape != np.shape(energy) or np.any(w < 0) or not np.isclose(w.sum(), 1):
            raise ValueError('Normalized nonnegative per-view law weights required')
    def score(t, nodes=64):
        j = gaussian_mixture_information(energy, t, nodes)
        a, b = float(wa @ j), float(wb @ j)
        return fraction * (1-t) * a - (1-fraction) * t * b, a, b
    root = brentq(lambda t: score(t)[0], 1e-7, 1-1e-7, xtol=1e-10)
    value, a, b = score(root)
    checked, aa, bb = score(root, 128)
    scale = fraction * (1-root) * aa + (1-fraction) * root * bb
    normalized = abs(checked) / scale if scale > 0 else float('inf')
    weak_a, weak_b = float(wa @ energy), float(wb @ energy)
    weak = fraction * weak_a / (fraction * weak_a + (1-fraction) * weak_b)
    return dict(fraction=root, bias=root-fraction, score_64=value, score_128=checked,
                normalized_score_check=normalized, quadrature_check_passed=bool(normalized <= 1e-6),
                information_a_64=a, information_b_64=b,
                information_a_128=aa, information_b_128=bb,
                information_odds_factor=a/b, weak_odds_factor=weak_a/weak_b,
                weak_fraction=weak, weak_bias=weak-fraction)


def histogram_law_weights(training_bins, metadata_counts):
    """Piecewise-Haar law evaluated by conditional within-bin sample means."""
    c = np.asarray(metadata_counts, float)
    bins = np.asarray(training_bins, int)
    n = np.bincount(bins, minlength=len(c))
    if np.any(c < 0) or c.sum() <= 0 or len(n) != len(c):
        raise ValueError('Invalid histogram')
    if np.any((c > 0) & (n == 0)):
        raise ValueError('Occupied metadata cell without a training rotation')
    p = c / c.sum()
    per_cell = np.divide(p, n, out=np.zeros_like(p), where=n > 0)
    return per_cell[bins]


def stratified_linear_summary(values, bins, metadata_counts):
    """MC SE conditional on cell counts and fixed metadata; not biological error."""
    x = np.asarray(values, float)
    bins = np.asarray(bins, int)
    c = np.asarray(metadata_counts, float)
    w = histogram_law_weights(bins, c)
    if not np.isfinite(x).all():
        return dict(mean=None, mc_standard_error=None, finite=False)
    n = np.bincount(bins, minlength=len(c))
    p = c / c.sum()
    sx = np.bincount(bins, weights=x, minlength=len(c))
    means = np.divide(sx, n, out=np.zeros_like(sx), where=n > 0)
    scatter = np.bincount(bins, weights=(x-means[bins])**2, minlength=len(c))
    if np.any((p > 0) & (n < 2)):
        se = None
    else:
        mean_var = np.divide(scatter, n*(n-1), out=np.zeros_like(scatter), where=n > 1)
        se = float(np.sqrt(np.sum(p*p*mean_var)))
    return dict(mean=float(w @ x), mc_standard_error=se, finite=True)
