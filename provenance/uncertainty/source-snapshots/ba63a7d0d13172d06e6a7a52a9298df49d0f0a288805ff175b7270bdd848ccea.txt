"""Adaptive reduced-space cubic design, with a separate final spectral audit.

Subgradient enrichment and cutting planes are classical optimization tools.
Approximate Ritz values guide selection; they are not upper certificates.
"""
import time
import numpy as np
from scipy.stats import norm
from .uq_cubic_subspace import ReducedCubicDesign, orthonormal_basis
from .uq_pose_optimization import gram_spectral_power
from .uq_random_spectral import gaussian_power_upper
from .uq_intervals import bias_aware_half_width_stable


def enrich_basis(basis, directions, tolerance=1e-10):
    """Preserve the existing span and append independent normalized directions."""
    q = np.asarray(basis, float).copy(); records = []
    for direction in directions:
        direction = np.asarray(direction, float)
        if direction.shape != (q.shape[0],) or not np.isfinite(direction).all():
            raise ValueError('Finite compatible enrichment direction required')
        length = float(np.linalg.norm(direction)); residual = direction.copy()
        for _ in range(2): residual -= q@(q.T@residual)
        remaining = float(np.linalg.norm(residual))
        added = length > 0 and remaining > tolerance*length and q.shape[1] < q.shape[0]
        records.append(dict(direction_norm=length, orthogonal_norm=remaining, added=bool(added)))
        if added: q = np.column_stack([q, residual/remaining])
    np.testing.assert_allclose(q.T@q, np.eye(q.shape[1]), atol=1e-9, rtol=1e-9)
    return q, records


def optimize_enriched_cubic(objective, basis, initial_weights, *, optimization_seed, certificate_seed,
                           outer_rounds=3, evaluations_per_round=6, design_seconds=7200,
                           separation_tolerance=1e-3, modes=3, ritz_tolerance=1e-4,
                           ritz_subspace=17, ritz_maxiter=150, power_probes=4, power_iterations=40,
                           alpha=.05/12, numerical_delta=1e-6/12, dense_oracle=False,
                           callback=None, checkpoint_callback=None):
    if optimization_seed == certificate_seed or min(optimization_seed, certificate_seed) < 0:
        raise ValueError('Distinct nonnegative design and certificate seeds required')
    if outer_rounds < 1 or evaluations_per_round < 2 or design_seconds <= 0:
        raise ValueError('Positive design budgets and at least two evaluations per round required')
    if not 0 < numerical_delta < alpha < .5 or abs(objective.z-norm.isf((alpha-numerical_delta)/2)) > 1e-12:
        raise ValueError('Inconsistent error budget')
    if objective.smoothing != 0: raise ValueError('Use an unsmoothed guide objective')
    begin = time.perf_counter(); q, _ = orthonormal_basis(basis)
    initial = np.asarray(initial_weights, float)
    np.testing.assert_allclose(q@(q.T@initial), initial, atol=1e-10, rtol=1e-9)
    residual, _, _ = objective.density(initial)
    ridge = float(objective.z*residual/(objective.B*np.linalg.norm(initial)))
    if not np.isfinite(ridge) or ridge <= 0: raise ValueError('Invalid fixed enrichment preconditioner')
    mode = np.random.default_rng(optimization_seed).normal(size=objective.op.shape[1]); mode /= np.linalg.norm(mode)
    best = dict(value=np.inf); history = []; rounds = []; stopped = 'round_budget'

    def announce(row):
        if callback: callback(row)

    # First evaluate in the inherited span, then enrich before each new master.
    for outer in range(outer_rounds+1):
        if outer and time.perf_counter()-begin >= design_seconds:
            stopped = 'design_wall_budget'; break
        enrichment = []
        if outer:
            gradient = best['gradient']
            directions = [-gradient, -gram_spectral_power(objective.gram, gradient, ridge, -1.)]
            q, enrichment = enrich_basis(q, directions)
        reduced = ReducedCubicDesign(objective, q)
        x = q.T@(initial if not best else best.get('weights', initial))
        cuts = []; master = None; round_history = []
        limit = 1 if outer == 0 else evaluations_per_round
        for local in range(limit):
            if history and time.perf_counter()-begin >= design_seconds:
                stopped = 'design_wall_budget'; break
            stamp = time.perf_counter()
            sigma, newcuts, spatial, mode, singular = reduced.spectral_directions(x, initial=mode, modes=modes,
                tolerance=ritz_tolerance, maxiter=ritz_maxiter, subspace=ritz_subspace, dense=dense_oracle)
            weights = q@x
            value, gradient, terms, support = objective.evaluate(weights, spatial)
            restricted_value, restricted_terms = reduced.value(x, sigma)
            row = dict(evaluation=len(history)+1, outer_round=outer, local_evaluation=local+1, rank=q.shape[1],
                guide_objective=value, restricted_guide_objective=restricted_value, terms=terms,
                restricted_terms=restricted_terms, ritz_singular_values=singular, master=master,
                coefficients=x.tolist(), seconds=time.perf_counter()-stamp, elapsed_seconds=time.perf_counter()-begin)
            history.append(row); round_history.append(row); cuts.extend(newcuts)
            # Cross-basis selection uses the original unsmoothed guide, not a
            # basis-dependent triangle majorant or unevaluated master proposal.
            if value < best['value']:
                best.update(value=value, weights=weights.copy(), gradient=gradient.copy(), support=support.copy(),
                            evaluation=row['evaluation'])
                if checkpoint_callback: checkpoint_callback(weights.copy(), row)
            announce(row)
            if outer == 0: break
            proposal, master = reduced.solve_master(cuts)
            row['next_master'] = master
            gap = max(0., restricted_value-master['objective'])/max(restricted_value, np.finfo(float).tiny)
            row['restricted_guide_gap'] = gap
            announce(dict(stage='enriched_master', outer_round=outer, local_evaluation=local+1, restricted_guide_gap=gap, **master))
            # Always evaluate at least one proposal in a newly enlarged span.
            if local >= 1 and gap <= separation_tolerance: break
            x = proposal
        rounds.append(dict(outer_round=outer, rank=q.shape[1], enrichment=enrichment, evaluations=len(round_history),
            diagnostics=reduced.diagnostics))
        if stopped == 'design_wall_budget': break
    if not np.isfinite(best['value']): raise RuntimeError('No feasible candidate evaluated')
    optimization_seconds = time.perf_counter()-begin
    objective.op.set_weights(best['weights'])
    announce(dict(stage='fresh_spectral_audit', evaluations=len(history), selected_evaluation=best['evaluation']))
    spectral = gaussian_power_upper(objective.op.spatial_gram, objective.op.shape[0], numerical_delta,
        power_probes, power_iterations, certificate_seed, callback=callback)
    audit = objective.audit(best['weights'], best['support'], spectral['eigenvalue_upper'])
    audit.update(weights=best['weights'], enriched_basis=q,
        half_width=bias_aware_half_width_stable(audit['noise_sd'], audit['bias'], alpha-numerical_delta),
        optimization_history=history, enrichment_rounds=rounds, selected_evaluation=best['evaluation'],
        selected_approximate_objective=best['value'], stopping_reason=stopped,
        optimization_seconds=optimization_seconds, certification_seconds=time.perf_counter()-begin-optimization_seconds,
        spectral_upper_certificate=spectral, preconditioner_ridge=ridge,
        optimization_seed=optimization_seed, certificate_seed=certificate_seed,
        alpha_total=alpha, alpha_noise=alpha-numerical_delta, alpha_numerical=numerical_delta,
        gap_scope='Restricted Ritz/conic gaps guide enrichment only. relative_sum_gap is the separately audited full-space triangle surrogate gap.',
        numerical_scope='Real-arithmetic integration inequalities; ordinary-floating-point and NUFFT guards, not validated arithmetic.')
    return audit
