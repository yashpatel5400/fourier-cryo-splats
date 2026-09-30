import importlib.util
from pathlib import Path


def test_relion_model_diagnostics_preserve_missing_and_unconverged_values(tmp_path):
    # A realistic scalar-plus-loop STAR file, independent of the parser's writer.
    star = '''data_model_general
_rlnCurrentResolution 96.48
_rlnAveragePmax 0.197123

data_model_classes
loop_
_rlnClassDistribution #1
_rlnAccuracyRotations #2
_rlnAccuracyTranslationsAngst #3
_rlnEstimatedResolution #4
1.0 6.66 7.680712 96.48
'''
    (tmp_path/'initial_it030_model.star').write_text(star)
    path = Path(__file__).resolve().parents[1]/'scripts/evaluate_relion_baseline.py'
    spec = importlib.util.spec_from_file_location('relion_evaluation', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    diagnostics = module.model_diagnostics(tmp_path)
    row = diagnostics['records'][0]
    assert row['general']['rlnCurrentResolution'] == 96.48
    assert row['general']['rlnPixelSize'] is None
    assert row['classes'][0]['rlnAccuracyRotations'] == 6.66
    assert row['classes'][0]['rlnAccuracyTranslationsAngst'] == 7.680712
    assert row['classes'][0]['rlnOverallFourierCompleteness'] is None
    assert diagnostics['units']['rlnAccuracyTranslationsAngst'] == 'angstrom'
    assert 'converged' not in row
    assert len(row['sha256']) == 64
