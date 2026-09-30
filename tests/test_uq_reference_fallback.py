from fourier_splats.uq_intervals import reference_interval_summary
import pytest


def test_no_data_choice_changes_center_and_variance_not_only_width():
    fallback=reference_interval_summary(2.,20.,.1,1.,4.,3.)
    assert fallback['uses_no_data']
    assert fallback['expected_center']==1. and fallback['selected_noise_sd']==0.
    assert fallback['analytic_coverage']==1. and fallback['correct_sign_probability']==0.
    # A refined bound can switch an originally no-data fit back to its raw
    # affine estimator; the original selected center must not be reused.
    refined=reference_interval_summary(2.,20.,.1,1.,2.,3.)
    assert not refined['uses_no_data']
    assert refined['expected_center']==20. and refined['selected_noise_sd']==.1
    assert refined['analytic_coverage']==0. and refined['correct_sign_probability']==1.
    with pytest.raises(ValueError,match='explicitly supplied pilot'):
        reference_interval_summary(2.,20.,.1,None,4.,3.)
    assert not reference_interval_summary(2.,20.,.1,None,2.,3.)['uses_no_data']
