"""Classical trust-region witnesses and residual-controlled dual upper bounds.

The probability, if any, belongs to a supplied spectral-upper event. These
routines create no new coverage guarantee and do not validate floating point.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def dense_quadratic_witness(gram, cross, radius):
    """Feasible guide for max v'Gv-2b'v on a ball, including the hard case.

    A secular-equation eigensolve is classical. Returned values are numerical
    lower diagnostics, never upper certificates. Feasibility is rechecked.
    """
    g = np.asarray(gram, float); b = np.asarray(cross, float); L = float(radius)
    if g.shape != (b.size, b.size) or b.ndim != 1 or b.size == 0 or not np.isfinite(g).all() or not np.isfinite(b).all() or not np.isfinite(L) or L < 0:
        raise ValueError('Finite symmetric Gram, compatible vector and nonnegative radius required')
    if not np.allclose(g, g.T, rtol=1e-10, atol=1e-12): raise ValueError('Symmetric Gram required')
    g = (g+g.T)/2
    if L == 0: return np.zeros_like(b), {'hard_case': False, 'value': 0., 'radius': L}
    d, u = np.linalg.eigh(g); top = float(d[-1]); c = u.T@b
    scale = max(float(np.linalg.norm(g, 2)), float(np.linalg.norm(b)/L), np.finfo(float).tiny)
    if d[0] < -1e-10*scale: raise ValueError('PSD Gram required')
    gaps = top-d; mask = gaps <= 64*np.finfo(float).eps*scale
    v = np.zeros_like(c); v[~mask] = -c[~mask]/gaps[~mask]
    hard = np.linalg.norm(c[mask]) <= 64*np.finfo(float).eps*scale*L and np.linalg.norm(v) <= L
    if hard:
        v[np.flatnonzero(mask)[-1]] = np.sqrt(max(0., L*L-float(v@v)))
        shift = 0.
    else:
        # Search in a dimensionless positive shift to avoid an absolute root
        # tolerance that is much larger than a small empirical Gram eigenvalue.
        def secular(delta):
            return float(np.linalg.norm(c/(gaps+scale*delta))/L-1.)
        low = np.finfo(float).tiny
        high = max(1., float(np.linalg.norm(c)/(L*scale)))
        while secular(high) > 0: high *= 2.
        if secular(low) <= 0:
            # Numerically near-hard: use a feasible direction; no upper claim.
            delta = low
        else:
            delta = brentq(secular, low, high, xtol=1e-14, rtol=1e-14)
        shift = float(scale*delta); v = -c/(gaps+shift)
    v = u@v
    length = float(np.linalg.norm(v))
    # For a PSD maximization the selected direction has b'v<=0, so
    # extending it to the boundary cannot reduce the objective. This also
    # avoids loss of near-hard-case accuracy from a tiny secular root error.
    if length > 0: v *= np.nextafter(L, 0.)/length
    return v, {'hard_case': bool(hard), 'value': float(v@g@v-2*b@v),
        'radius': L, 'witness_norm': float(np.linalg.norm(v)), 'shift_above_ritz': shift,
        'largest_ritz_value': top, 'scope': 'Feasible numerical guide, not an upper certificate.'}


def krylov_cache(matvec, start, steps=30, extras=()):
    """Fully reorthogonalized operator cache; store actual GQ, not a recurrence.

    Additional starts help the guide see a leading eigenspace orthogonal to b.
    Every retained column receives one explicitly evaluated operator action.
    """
    b = np.asarray(start, float)
    if b.ndim != 1 or not b.size or not np.isfinite(b).all() or steps < 1:
        raise ValueError('Finite nonempty start and positive step count required')
    vectors = []; images = []; pending = [b.copy(), *[np.asarray(x, float).copy() for x in extras]]
    while pending and len(vectors) < min(steps, b.size):
        raw = pending.pop(0)
        if raw.shape != b.shape or not np.isfinite(raw).all(): raise ValueError('Compatible finite extra starts required')
        norm = float(np.linalg.norm(raw)); q = raw.copy()
        if vectors:
            basis = np.column_stack(vectors)
            for _ in range(2): q -= basis@(basis.T@q)
        length = float(np.linalg.norm(q))
        if norm == 0 or length <= 1e-12*norm: continue
        q /= length; gq = np.asarray(matvec(q), float)
        if gq.shape != b.shape or not np.isfinite(gq).all(): raise FloatingPointError('Invalid Gram action')
        vectors.append(q); images.append(gq); pending.append(gq.copy())
    q = np.column_stack(vectors) if vectors else np.zeros((b.size, 0))
    gq = np.column_stack(images) if images else np.zeros_like(q)
    reduced = q.T@gq
    return q, gq, {'rank': q.shape[1], 'requested_steps': steps,
        'orthogonality_defect': float(np.linalg.norm(q.T@q-np.eye(q.shape[1]))),
        'reduced_symmetry_defect': float(np.linalg.norm(reduced-reduced.T))}


def resolvent_upper(cross, spectral_upper, shift, basis, images):
    """Bound b'(lambda I-G)^-1 b from an explicit Galerkin residual.

    Assumes G symmetric PSD and lambda_max(G)<=spectral_upper. Let A=lambda I-G,
    x=Q y and r=b-Ax. The exact identity is b'A^-1 b=2b'x-x'Ax+r'A^-1 r.
    The last term is at most ||r||²/(lambda-spectral_upper). The identity
    holds for any x, so inaccurate reduced solves cannot invalidate it in
    real arithmetic. Floating-point roundoff is not enclosed here.
    """
    b = np.asarray(cross, float); q = np.asarray(basis, float); gq = np.asarray(images, float)
    U = float(spectral_upper); lam = float(shift)
    if b.ndim != 1 or q.ndim != 2 or q.shape[0] != b.size or q.shape != gq.shape or not all(np.isfinite(x).all() for x in [b,q,gq]) or not np.isfinite([U,lam]).all() or U < 0 or lam <= U:
        raise ValueError('Finite compatible cache and lambda > U >= 0 required')
    if q.shape[1]:
        reduced = lam*(q.T@q)-.5*(q.T@gq+gq.T@q)
        y = np.linalg.solve(reduced, q.T@b)
        x = q@y; ax = lam*x-gq@y
    else:
        x = np.zeros_like(b); ax = x.copy()
    residual = b-ax
    base = float(2*b@x-x@ax); correction = float(residual@residual/(lam-U))
    value = base+correction
    if not np.isfinite(value) or value < -1e-10*max(1.,abs(base),correction):
        raise FloatingPointError('Invalid resolvent upper value')
    return {'upper': float(max(0.,value)), 'galerkin_term': base, 'residual_correction': correction,
        'residual_norm': float(np.linalg.norm(residual)), 'lambda': lam, 'spectral_upper': U}


def quadratic_dual_upper(cross, radius, spectral_upper, basis, images):
    """Minimize valid dual upper values; all use the same spectral event.

    A fixed log grid and bounded scalar refinement only select among valid
    upper bounds. Neither solver convergence nor the Krylov rank is a premise.
    """
    b = np.asarray(cross, float); L = float(radius); U = float(spectral_upper)
    if not np.isfinite(L) or L <= 0 or not np.isfinite(U) or U < 0: raise ValueError('Positive radius and nonnegative upper required')
    scale = max(U, float(np.linalg.norm(b)/L), 1e-300)
    records = []
    def evaluate(logshift):
        lam = U+scale*np.exp(logshift)
        item = resolvent_upper(b,U,lam,basis,images)
        item['quadratic_upper'] = float(lam*L*L+item['upper'])
        item['log_shift'] = float(logshift); records.append(item)
        return item['quadratic_upper']
    grid = np.linspace(-20., 10., 61)
    values = np.array([evaluate(x) for x in grid]); j = int(np.argmin(values))
    opt = minimize_scalar(evaluate, bounds=(grid[max(0,j-1)],grid[min(len(grid)-1,j+1)]),
        method='bounded', options={'xatol': 1e-7, 'maxiter': 100})
    selected = min(records, key=lambda row: row['quadratic_upper'])
    triangle_quadratic = float(U*L*L+2*L*np.linalg.norm(b))
    return {'quadratic_upper': min(selected['quadratic_upper'],triangle_quadratic),
        'resolvent_quadratic_upper': selected['quadratic_upper'], 'old_cross_quadratic_upper': triangle_quadratic,
        'selected': selected, 'evaluations': records, 'scalar_success': bool(opt.success),
        'scope': 'Conditional on the supplied spectral event, in real arithmetic; no validated floating-point enclosure.'}
