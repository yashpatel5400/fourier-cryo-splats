"""Explicit centered-cube support assumptions and an outside-support allowance.

These change the density class. They are not estimates of molecular support.
The coordinate change is an L2 isometry, and the outside allowance is a
conservative continuous Fourier bound, not a finite-voxel guarantee.
"""
import numpy as np


def rescale_supported_problem(k, q, ctf, centers, signs, width, shift, side):
    """Map a density on [-side/2, side/2]^3 to the unit cube isometrically.

    rho_tilde(u)=side**(3/2) rho(side*u), so the observation transfer gains
    side**(3/2). A normalized Gaussian target gains side**(-3/2).
    Translation remains physical: (side*q)*(shift/side)=q*shift.
    """
    side = float(side)
    if not np.isfinite(side) or not 0 < side <= 1:
        raise ValueError('Support side must be in (0, 1]')
    if width <= 0 or shift < 0:
        raise ValueError('Positive target width and nonnegative shift required')
    return {'k': np.asarray(k)*side, 'q': np.asarray(q)*side,
            'ctf': np.asarray(ctf)*side**1.5,
            'centers': np.asarray(centers)/side,
            'signs': np.asarray(signs)/side**1.5,
            'width': float(width/side), 'shift': float(shift/side),
            'side': side}


def restrict_cell_density(coefficients, box, side):
    """Exact restriction of orthonormal cell coefficients to aligned cube edges.

    No interpolation, renormalization, clipping, or positivity projection is
    performed. Returned coefficients represent the isometrically rescaled
    density; their Euclidean norm is its continuous L2 norm.
    """
    count = int(round(box*side))
    if not 0 < side <= 1 or abs(count-box*side) > 1e-10 or (box-count) % 2:
        raise ValueError('Support boundaries must coincide with centered cell edges')
    start = (box-count)//2
    a = np.asarray(coefficients).reshape(box, box, box)
    cropped = a[start:start+count, start:start+count, start:start+count].copy()
    tail = np.sqrt(max(0., float(np.sum(a*a)-np.sum(cropped*cropped))))
    return cropped.ravel(), count, float(tail)


def outside_tail_bias_upper(k, q, ctf, weights, noise, angle, shift,
                            nominal_residual_norm, tail_radius):
    """Bound an arbitrary outside density with supplied continuous L2 radius.

    For each Fourier column, symmetry of the cube gives
    ||exp(i phi_u)-exp(i phi_0)||_2 <= min(2, 2*pi*sqrt(
        angle**2*||k||**2/12 + shift**2*||q||**2)).
    The rotation chord is at most angle*||k||. The triangle inequality then
    bounds the full perturbed residual, hence also its restriction outside.
    No assumption that the tail is spatially white or independent is made.
    """
    if noise <= 0 or min(angle, shift, nominal_residual_norm, tail_radius) < 0:
        raise ValueError('Invalid noise scale, radii or residual norm')
    k, q, ctf = map(np.asarray, (k, q, ctf))
    n, nq = ctf.shape
    w = np.asarray(weights).reshape(n, 2*nq)
    amplitude = np.hypot(w[:, :nq], w[:, nq:])*abs(ctf)/noise
    column = np.minimum(2., 2*np.pi*np.sqrt(
        angle**2*np.sum(k*k, axis=-1)/12 + shift**2*np.sum(q*q, axis=-1)))
    pose = float(np.sum(amplitude*column))
    return {'bias_upper': float(tail_radius*(nominal_residual_norm+pose)),
            'nominal_residual_upper': float(nominal_residual_norm),
            'pose_field_upper': pose, 'tail_radius_supplied': float(tail_radius)}
