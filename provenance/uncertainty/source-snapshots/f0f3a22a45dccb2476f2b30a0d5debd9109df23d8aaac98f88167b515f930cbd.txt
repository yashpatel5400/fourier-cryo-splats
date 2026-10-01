"""Aggregate a paired study without dropping failed or unfinished replicates."""
from itertools import product
import numpy as np
from .uq_end_to_end import binomial_interval


TEMPLATES = ('oracle_reference', 'independent_pilot')
TARGETS = ('center', 'contrast')
METHODS = ('fixed_folded', 'fixed_sum', 'deterministic_pose', 'mixed_common0',
           'mixed_common01', 'gaussian_tau2_fixed', 'gaussian_tau2_pose',
           'gaussian_tau1_fixed', 'gaussian_tau1_pose')
IMAGES = ('same_image', 'independent_image')
KEYS = ('template', 'target', 'method', 'image_mode')


def distribution(values):
    a = np.asarray([v for v in values if v is not None], float)
    if not len(a):
        return dict(count=0)
    if not np.isfinite(a).all():
        raise ValueError('Nonfinite summary input')
    return dict(count=len(a), mean=float(a.mean()), minimum=float(a.min()),
                median=float(np.median(a)), p05=float(np.quantile(a, .05)),
                p95=float(np.quantile(a, .95)), maximum=float(a.max()))


def aggregate_replicates(records, planned=200):
    """Intervals within a replicate stay paired; each method has planned trials.

    Completed but nonconverged numerical procedures remain in the evaluation.
    An exception/incomplete replicate is an algorithmic failure, not a reason
    to remove a trial. Planned-denominator coverage counts its result as a
    failure. Exact binomial intervals are emitted only when all trials were
    attempted. For unfinished studies, bounds show all missing outcomes.
    """
    if not isinstance(planned, int) or planned < 1:
        raise ValueError('Positive planned replication count required')
    ids = [r['replicate'] for r in records]
    if len(set(ids)) != len(ids) or any(not isinstance(i, int) or not 0 <= i < planned for i in ids):
        raise ValueError('Duplicate or out-of-range replicate')
    attempted = [r for r in records if r.get('complete') or r.get('error')]
    finished = len(attempted) == planned
    cells = {k: [] for k in product(TEMPLATES, TARGETS, METHODS, IMAGES)}
    for record in records:
        if not record.get('complete'):
            continue
        seen = set()
        for row in record['intervals']:
            key = tuple(row[k] for k in KEYS)
            if key not in cells or key in seen:
                raise ValueError('Unexpected or duplicated interval cell')
            seen.add(key)
            cells[key].append((record['replicate'], row))
        if seen != set(cells):
            raise ValueError('Completed replicate is missing declared interval cells')
    output = dict(planned_replicates=planned, attempted_replicates=len(attempted),
        completed_replicates=sum(bool(r.get('complete')) for r in records),
        failed_replicates=sum(bool(r.get('error')) and not r.get('complete') for r in records),
        unattempted_replicates=planned-len(attempted), complete=finished, groups=[], paired_image_controls=[],
        scope='One Bernoulli trial per independent simulated dataset and procedure. Images, templates, targets and methods are paired. Width summaries use available completed replicates only; failed trials remain in the planned coverage denominator.')
    boolean_fields = ('covered', 'raw_covered', 'correct_sign_exclusion', 'raw_correct_sign_exclusion', 'fallback')
    for key, entries in cells.items():
        rows = [r for _, r in entries]
        group = dict(zip(KEYS, key)); group['available_intervals'] = len(rows)
        group['missing_or_failed_intervals'] = planned-len(rows)
        for field in boolean_fields:
            successes = sum(bool(r[field]) for r in rows)
            group[field] = dict(successes=successes, planned_denominator=planned,
                fraction=successes/planned,
                exact_binomial_95=binomial_interval(successes, planned) if finished else None,
                unresolved_fraction_range=[successes/planned, (successes+planned-len(rows))/planned])
        group['relative_width_below_half_count'] = sum(r['relative_half_width'] < .5 for r in rows)
        for field in ('raw_half_width', 'half_width', 'raw_relative_half_width', 'relative_half_width',
                      'half_width_over_abs_center', 'raw_center', 'center'):
            group[field] = distribution([r[field] for r in rows])
        group['signed_error'] = distribution([r['center']-r['true_target'] for r in rows])
        group['raw_signed_error'] = distribution([r['raw_center']-r['true_target'] for r in rows])
        group['truth'] = distribution([r['true_target'] for r in rows])
        output['groups'].append(group)
    for template, target, method in product(TEMPLATES, TARGETS, METHODS):
        same = dict(cells[(template, target, method, 'same_image')])
        other = dict(cells[(template, target, method, 'independent_image')])
        if same.keys() != other.keys():
            raise ValueError('Paired image controls have inconsistent replication')
        table = {'both_cover': 0, 'same_only': 0, 'independent_only': 0, 'neither_covers': 0}
        for rep, a in same.items():
            b = other[rep]
            label = ('both_cover' if b['covered'] else 'same_only') if a['covered'] else (
                'independent_only' if b['covered'] else 'neither_covers')
            table[label] += 1
        output['paired_image_controls'].append(dict(template=template, target=target, method=method,
            available_pairs=len(same), **table,
            independent_minus_same_coverage=(table['independent_only']-table['same_only'])/planned,
            paired_center_difference=distribution([other[i]['raw_center']-same[i]['raw_center'] for i in same])))
    return output
