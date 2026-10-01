# Fresh constrained Fisher comparison: complete results

All three candidates use 65,536 new training views, 32,768 new calibration views with two groups of 32 noises, and 131,072 independent Haar test views. Matched and constrained Fisher scores share the same training mean gap; Fisher uses covariance shrinkage 0.1 without tuning. All twelve scores are retained. Seven calibration methods, five viewing caps, three particle counts, five deletions and three amplitudes are in the CSVs.

## Detection at n=10,000, amplitude one

Entries are matched / Fisher rejection projections. This table includes both paired variance and split CVaR; neither score nor calibration method is selected by an unadjusted minimum. Other methods, amplitudes, caps and particle counts remain in the full archive.

| Stack | Family | Bound | kappa | 10% deletion | 25% | 50% | 100% |
|---|---|---|---:|---:|---:|---:|---:|
| 10028 | power | view_variance | 1 | 1.48e-08 / 2.89e-08 | 6.58e-05 / 0.000156 | 0.216 / 0.345 | 1 / 1 |
| 10028 | power | view_variance | 1.1 | 1.6e-56 / 7.73e-53 | 2.48e-45 / 1.82e-41 | 1.08e-28 / 6.73e-25 | 1.96e-06 / 5.29e-05 |
| 10028 | power | cvar_split | 1 | 8.3e-09 / 1.46e-08 | 4.37e-05 / 9.74e-05 | 0.188 / 0.302 | 1 / 1 |
| 10028 | power | cvar_split | 1.1 | 8.29e-29 / 1.18e-27 | 4.14e-21 / 9.71e-20 | 1.29e-10 / 3.08e-09 | 0.55 / 0.712 |
| 10028 | power_bispectrum | view_variance | 1 | 6.64e-10 / 4.81e-09 | 2.61e-06 / 7.65e-05 | 0.031 / 0.269 | 1 / 1 |
| 10028 | power_bispectrum | view_variance | 1.1 | 1.14e-57 / 1.46e-55 | 1.06e-47 / 4.31e-43 | 2.42e-32 / 2.53e-26 | 5e-10 / 0.000102 |
| 10028 | power_bispectrum | cvar_split | 1 | 5.16e-10 / 2.99e-09 | 2.15e-06 / 5.53e-05 | 0.0283 / 0.243 | 1 / 1 |
| 10028 | power_bispectrum | cvar_split | 1.1 | 1.67e-31 / 1.13e-29 | 2.43e-24 / 6.8e-21 | 5.5e-14 / 4.2e-10 | 0.0407 / 0.762 |
| 10049 | power | view_variance | 1 | 0.000219 / 0.0003 | 0.0456 / 0.0603 | 0.909 / 0.905 | 1 / 1 |
| 10049 | power | view_variance | 1.1 | 1.85e-16 / 7.18e-16 | 1.29e-10 / 5.14e-10 | 0.000475 / 0.000591 | 0.978 / 0.987 |
| 10049 | power | cvar_split | 1 | 0.000101 / 0.000141 | 0.0294 / 0.0398 | 0.871 / 0.867 | 1 / 1 |
| 10049 | power | cvar_split | 1.1 | 5.1e-12 / 1.39e-11 | 3.22e-07 / 8.75e-07 | 0.0251 / 0.0276 | 1 / 1 |
| 10049 | power_bispectrum | view_variance | 1 | 5.71e-05 / 9.64e-05 | 0.0289 / 0.0513 | 0.899 / 0.928 | 1 / 1 |
| 10049 | power_bispectrum | view_variance | 1.1 | 1.02e-17 / 7.32e-17 | 3.26e-11 / 3.51e-10 | 0.000385 / 0.00104 | 0.97 / 0.983 |
| 10049 | power_bispectrum | cvar_split | 1 | 2.27e-05 / 3.61e-05 | 0.0171 / 0.0305 | 0.854 / 0.888 | 1 / 1 |
| 10049 | power_bispectrum | cvar_split | 1.1 | 5.17e-13 / 1.53e-12 | 1.2e-07 / 5.29e-07 | 0.0229 / 0.0365 | 0.999 / 1 |
| 10076 | power | view_variance | 1 | 5.34e-05 / 7.88e-05 | 0.000324 / 0.000454 | 0.00613 / 0.00565 | 0.155 / 0.159 |
| 10076 | power | view_variance | 1.1 | 3.22e-10 / 5.9e-10 | 5.46e-09 / 9.42e-09 | 7.57e-07 / 6.55e-07 | 0.000448 / 0.000479 |
| 10076 | power | cvar_split | 1 | 1.93e-05 / 2.68e-05 | 0.00013 / 0.000173 | 0.00303 / 0.0026 | 0.105 / 0.104 |
| 10076 | power | cvar_split | 1.1 | 3.51e-11 / 5.82e-11 | 6.95e-10 / 1.09e-09 | 1.3e-07 / 1e-07 | 0.000125 / 0.000124 |
| 10076 | power_bispectrum | view_variance | 1 | 4.61e-05 / 3.39e-05 | 0.000395 / 0.000259 | 0.00632 / 0.00438 | 0.18 / 0.17 |
| 10076 | power_bispectrum | view_variance | 1.1 | 2.56e-10 / 1.6e-10 | 7.49e-09 / 3.82e-09 | 7.95e-07 / 4.19e-07 | 0.000639 / 0.000556 |
| 10076 | power_bispectrum | cvar_split | 1 | 1.81e-05 / 1.2e-05 | 0.000174 / 0.000103 | 0.00332 / 0.00211 | 0.128 / 0.116 |
| 10076 | power_bispectrum | cvar_split | 1.1 | 3.15e-11 / 1.69e-11 | 1.1e-09 / 4.77e-10 | 1.52e-07 / 6.92e-08 | 0.000198 / 0.000158 |

## Pointwise difference intervals

Each row below compares 25% deletion at n=100,000, amplitude one, kappa=1. Intervals combine two exact 97.5% single-image intervals by a union bound and monotone binomial transforms. They have pointwise at least 95% coverage conditional on frozen training/calibration under the simulator. They do not quantify training/calibration variability and are not simultaneous. All 9,450 pairs and their joint event tables are retained.

| Stack | Family | Method | Matched | Fisher | Difference | Interval |
|---|---|---|---:|---:|---:|---|
| 10028 | power | view_variance | 1.34e-17 | 2.86e-15 | 2.84e-15 | [-4.01e-11, 2.45e-09] |
| 10028 | power | cvar_split | 5.68e-19 | 1.31e-16 | 1.3e-16 | [-3.41e-12, 2.33e-10] |
| 10028 | power_bispectrum | view_variance | 1.55e-27 | 2.69e-17 | 2.69e-17 | [-4.33e-19, 6.87e-11] |
| 10028 | power_bispectrum | cvar_split | 3.34e-28 | 2.86e-18 | 2.86e-18 | [-1.22e-19, 1.21e-11] |
| 10049 | power | view_variance | 0.042 | 0.0941 | 0.0521 | [-0.591, 0.74] |
| 10049 | power | cvar_split | 0.00821 | 0.026 | 0.0178 | [-0.33, 0.507] |
| 10049 | power_bispectrum | view_variance | 0.00844 | 0.0599 | 0.0515 | [-0.334, 0.657] |
| 10049 | power_bispectrum | cvar_split | 0.00095 | 0.0106 | 0.0097 | [-0.126, 0.366] |
| 10076 | power | view_variance | 3.03e-13 | 2.63e-12 | 2.32e-12 | [-8.05e-08, 3.93e-07] |
| 10076 | power | cvar_split | 1.12e-15 | 7.73e-15 | 6.61e-15 | [-1.2e-09, 5.2e-09] |
| 10076 | power_bispectrum | view_variance | 1.22e-12 | 7.02e-14 | -1.15e-12 | [-2.24e-07, 2.72e-08] |
| 10076 | power_bispectrum | cvar_split | 7.67e-15 | 1.83e-16 | -7.49e-15 | [-5.17e-09, 3.01e-10] |

## Correct-null controls and verification

Largest n=10,000 correct-null projection: 2.7347584e-05, interval [2.3639861e-06, 0.00023906012], at 10076, power_fisher, grouped_ratio, kappa=1.0, amplitude=1.1.
Largest correct-null actual group count: 2/128, pointwise interval [0.0018979, 0.055303]. Calibration is fixed, cells share draws, and the 1,260 null cells cannot be pooled as independent repetitions.

Independent replay checks 3 full feature covariances, 12 directions via a null-space optimization, 420 probability bounds, 1680 minimal critical values, 18900 projections and 6300 group vectors. Physical first views, first calibration cubics and all training noise draws are also checked. This does not certify all physical operator calls or rounding.

Whole-candidate testing still depends on a hand-specified deletion family and known transfer/noise/viewing assumptions. It does not localize a rejection to the nominated region or prove experimental regional-density coverage. This classical baseline does not supply a new optimization theorem or an ICML acceptance assessment.
