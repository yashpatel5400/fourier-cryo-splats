"""Acquisition-only fixed-grid diagnostics; no claim to identify pure noise."""
import numpy as np


def correlation(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Matching finite fields required')
    a, b = a-a.mean(), b-b.mean()
    denominator = np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.sum(a*b)/denominator) if denominator else None


def radial_power(field):
    """Unwindowed demeaned periodogram; shell energies sum to pixel variance."""
    x = np.asarray(field, float)
    if x.ndim != 2 or x.shape[0] != x.shape[1] or len(x) % 2 or not np.isfinite(x).all():
        raise ValueError('Finite square even-sized patch required')
    n = len(x); x = x-x.mean()
    f = np.fft.rfft2(x)
    power = abs(f)**2/n**4
    multiplicity = np.full(power.shape, 2.); multiplicity[:, [0, -1]] = 1.
    ky, kx = np.fft.fftfreq(n)*n, np.fft.rfftfreq(n)*n
    shell = np.floor(np.hypot(ky[:, None], kx[None, :])).astype(int)
    energy = np.bincount(shell.ravel(), weights=(power*multiplicity).ravel())
    modes = np.bincount(shell.ravel(), weights=multiplicity.ravel())
    return dict(shell_energy=energy, modes=modes,
                mean_power=np.divide(energy, modes, out=np.zeros_like(energy), where=modes > 0),
                variance=float(np.mean(x*x)))


def block_average(x, factor):
    x = np.asarray(x)
    if x.ndim != 2 or not isinstance(factor, int) or factor < 1 or any(n % factor for n in x.shape):
        raise ValueError('Block factor must divide both image dimensions')
    return x.reshape(x.shape[0]//factor, factor, x.shape[1]//factor, factor).mean(axis=(1, 3))
