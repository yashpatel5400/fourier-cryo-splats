"""Check that review compaction preserves unfavorable case records and identity."""
import importlib.util
from pathlib import Path

PATH = Path(__file__).resolve().parents[1]/'scripts/review_uq_candidate.py'
spec = importlib.util.spec_from_file_location('scientific_review_packet', PATH)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


def test_case_records_and_failures_are_not_pruned():
    cases = [{'complete': True, 'analytic_coverage': i/200, 'failure': i < 190} for i in range(200)]
    projected = module.compact_evidence({'records': cases})['records']
    restored = [dict(projected['common_fields'], **dict(zip(projected['columns'], row))) for row in projected['rows']]
    assert restored == cases
    assert projected['original_count'] == 200


def test_long_vectors_have_replayable_identity_and_visible_omission():
    a = list(range(1000)); result = module.compact_evidence(a)
    assert result['original_count'] == 1000
    assert result['minimum'] == 0 and result['maximum'] == 999
    assert result['first'] == [0, 1, 2] and result['last'] == [997, 998, 999]
    assert 'summarized' in result['review_projection']
    b = a.copy(); b[400] += 1
    assert result['canonical_json_sha256'] != module.compact_evidence(b)['canonical_json_sha256']


def test_short_outcome_vectors_and_boolean_counts_remain_clear():
    assert module.compact_evidence([.94, .99, .1]) == [.94, .99, .1]
    result = module.compact_evidence([True]*75+[False]*80)
    assert result['original_count'] == 155 and result['true_count'] == 75


def test_float_rounding_does_not_change_failure_flags_or_small_exponents():
    row = module.compact_evidence({'radius': 1.00000000001, 'in_class': False, 'error': 1.234567890123e-42, 'seed': 20260930})
    assert row['in_class'] is False and row['seed'] == 20260930
    assert abs(row['error']/1.234567890123e-42-1) < 1e-7


def test_flattened_baseline_deduplication_checks_every_original_row():
    import copy
    import json
    path = module.ROOT/'results/uncertainty/development/pilot-selected-fourier-summary-v1/summary.json'
    if not path.exists():
        import pytest
        pytest.skip('Completed baseline summary needed for integration check')
    raw = json.loads(path.read_text())
    projected = module.remove_duplicated_baseline_rows(path, raw)
    assert projected['records']['original_count'] == 960
    changed = copy.deepcopy(raw); changed['records'][123]['analytic_fixed_signal_coverage'] = -1.
    import pytest
    with pytest.raises(ValueError, match='not a duplicate'):
        module.remove_duplicated_baseline_rows(path, changed)


def test_read_only_selection_uses_file_role_not_favorable_outcomes():
    assert module.inline_with_read_tools('results/uncertainty/development/study/summary.json')
    assert module.inline_with_read_tools('results/uncertainty/development/noise-scale-calibration.json')
    assert module.inline_with_read_tools('paper/main.tex')
    assert not module.inline_with_read_tools('src/fourier_splats/uq_continuous.py')
    assert not module.inline_with_read_tools('scripts/apply_uq_fresh_noise.py')
    assert not module.inline_with_read_tools('tests/test_uq_continuous.py')
    assert not module.inline_with_read_tools('results/uncertainty/development/study/success.json')
    assert not module.inline_with_read_tools('results/uncertainty/development/study/failure.json')


def test_snapshot_copy_accepts_only_exact_content_addressed_paths():
    good = 'provenance/uncertainty/source-snapshots/'+'a'*64+'.txt'
    bad = 'provenance/uncertainty/source-snapshots/../../'+'a'*64+'.txt'
    assert list(module.snapshot_references({'sources': [{'snapshot': good}, {'snapshot': bad}]})) == [good]
