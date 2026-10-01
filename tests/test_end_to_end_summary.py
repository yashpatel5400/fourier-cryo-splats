from itertools import product
import pytest
from fourier_splats.uq_end_to_end import observed_interval
from fourier_splats.uq_end_to_end_summary import (
    aggregate_replicates, TEMPLATES, TARGETS, METHODS, IMAGES)


def record(i, same_covers=True):
    rows = []
    for template, target, method, image in product(TEMPLATES, TARGETS, METHODS, IMAGES):
        center = 1. if image == 'independent_image' or same_covers else 3.
        rows.append(dict(template=template, target=target, method=method, image_mode=image,
            **observed_interval(center, .5, 0., 2., 1., False)))
    return dict(replicate=i, complete=True, intervals=rows)


def test_failed_dataset_stays_in_planned_denominator():
    result = aggregate_replicates([record(0), dict(replicate=1, complete=False, error='solver failure')], 2)
    assert result['complete'] and result['failed_replicates'] == 1
    assert len(result['groups']) == 72
    for group in result['groups']:
        assert group['covered']['successes'] == 1
        assert group['covered']['fraction'] == .5
        assert group['covered']['exact_binomial_95'][1] > .9
        assert group['half_width']['count'] == 1


def test_partial_study_has_no_completed_monte_carlo_interval():
    result = aggregate_replicates([record(0)], 2)
    assert not result['complete']
    assert result['groups'][0]['covered']['exact_binomial_95'] is None
    assert result['groups'][0]['covered']['unresolved_fraction_range'] == [.5, 1.]


def test_paired_discordance_and_invalid_replications():
    result = aggregate_replicates([record(0, False), record(1)], 2)
    assert all(p['independent_only'] == 1 and p['both_cover'] == 1 for p in result['paired_image_controls'])
    assert all(p['independent_minus_same_coverage'] == .5 for p in result['paired_image_controls'])
    with pytest.raises(ValueError, match='Duplicate'):
        aggregate_replicates([record(0), record(0)], 2)
    missing = record(0); missing['intervals'].pop()
    with pytest.raises(ValueError, match='missing'):
        aggregate_replicates([missing], 1)
