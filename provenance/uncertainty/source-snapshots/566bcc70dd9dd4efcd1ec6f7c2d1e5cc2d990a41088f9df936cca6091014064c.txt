"""Higher-order remainder diagnostics; not a higher-order confidence interval.

The omitted lower-order pose fields must also be retained before using a new
remainder in an interval. Floating-point guards are not interval arithmetic.
"""
import math
import numpy as np
from .uq_ball_remainder import particle_ball_remainder, ball_sobolev_norms


def partial_bell(derivative_bounds):
    """Nonnegative partial exponential Bell polynomials, including B[0,0]=1."""
    values = np.asarray(derivative_bounds, float)
    if values.ndim < 1 or not len(values) or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Finite nonnegative derivative bounds required')
    order = len(values)
    table = np.zeros((order+1, order+1)+values.shape[1:])
    table[0, 0] = 1
    for n in range(1, order+1):
        for j in range(1, n+1):
            for i in range(1, n-j+2):
                table[n, j] += math.comb(n-1, i-1)*values[i-1]*table[n-i, j-1]
    if not np.isfinite(table).all():
        raise FloatingPointError('Nonfinite Bell polynomial')
    return table


def higher_order_particle_remainders(k, q, coefficients, angle, shift,
                                      degrees=(2, 3, 4, 5), joint_ball=True, domain='cube'):
    degrees = tuple(degrees)
    if not degrees or any(not isinstance(d, (int, np.integer)) or not 2 <= d <= 7 for d in degrees):
        raise ValueError('Taylor degrees must be integers from two through seven')
    if len(set(degrees)) != len(degrees):
        raise ValueError('Repeated Taylor degree')
    # Reuse the checked geometry and embedding bounds, but not its cubic term.
    base = particle_ball_remainder(k, q, coefficients, angle, shift,
                                   joint_ball=joint_ball, domain=domain)
    k = np.asarray(k, float); q = np.asarray(q, float); c = np.asarray(coefficients, complex)
    order = max(degrees)+1
    sobolev, diagnostics = ball_sobolev_norms(k, c, base['domain_radius'],
                                            orders=tuple(range(1, order+1)), domain=domain)
    path = np.array([base['path_speed_bound']]+[angle**j*np.sqrt(3)/2 for j in range(2, order+1)])
    bell = partial_bell(path)
    embedding = np.linalg.lstsq(q, k, rcond=None)[0]
    eps = 2*np.pi*shift*np.linalg.norm(q-k@np.linalg.pinv(embedding), axis=1)
    phase = path[:, None]*(2*np.pi*np.linalg.norm(k, axis=1))[None]
    complete_phase = partial_bell(phase).sum(axis=1)
    records = []
    for degree in degrees:
        n = degree+1
        main = float(bell[n, 1:n+1]@sobolev[:n]/math.factorial(n))
        residual_per_frequency = eps*complete_phase[n]
        for j in range(1, n+1):
            residual_per_frequency += math.comb(n, j)*eps**j*complete_phase[n-j]
        residual = float(abs(c)@residual_per_frequency/math.factorial(n))
        if not np.isfinite(main+residual):
            raise FloatingPointError('Nonfinite higher-order remainder')
        records.append({'taylor_degree': int(degree), 'main_remainder': main,
                        'embedding_residual_remainder': residual,
                        'field_remainder': main+residual})
    degree_two = next((r['field_remainder'] for r in records if r['taylor_degree'] == 2), None)
    return {'records': records, 'domain': domain, 'domain_radius': base['domain_radius'],
            'path_derivative_bounds': path.tolist(), 'sobolev_norms': sobolev.tolist(),
            'maximum_phase_residual': float(np.max(eps)),
            'degree_two_existing_bound': base['field_remainder'],
            'degree_two_absolute_difference': None if degree_two is None else abs(degree_two-base['field_remainder']),
            'maximum_relative_sobolev_pad': max(d['roundoff_pad']/max(abs(d['squared_unpadded']), np.finfo(float).tiny) for d in diagnostics),
            'scope': 'Remainder alone. Higher pose-polynomial terms have not been included in a confidence interval.'}
