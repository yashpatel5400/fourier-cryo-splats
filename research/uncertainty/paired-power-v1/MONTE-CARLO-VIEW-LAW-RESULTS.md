# Monte Carlo candidate validation: bounded-view feasibility results

The frozen study completes in 368.846 seconds. Each nonzero stack uses 131,072 independent continuous Haar calibration draws and 131,072 held-out draws. The true and region-removed maps share rotations/noise within a stage. Calibration and held-out stages are independent. Four zero 10076 directions are retained as uninformative, with no projections or random draws. The study has no new experimental images or reconstruction results.

The test assumes a known proper Gaussian noise/CTF model, an amplitude interval and a specified bound kappa on the viewing density relative to Haar measure. Its alpha=.05 guarantee includes delta=.001 calibration failure; it does not condition on arbitrary realized poses or establish those assumptions from experimental data. Each candidate/test is assessed separately.

## All nonzero contrasts at 100,000 tested particles

Entries are projected rejection probabilities on the true-map Haar alternative, using the independently calibrated region-removed candidate. Brackets are pointwise 95% intervals obtained from held-out single-particle probabilities. They are not observed repeated-dataset power or simultaneous intervals.

| Stack | Contrast | kappa=1 | kappa=1.01 | kappa=1.1 | kappa=2 |
|---|---|---:|---:|---:|---:|
| 10028 | power_fixed | 1 [1, 1] | 1 [1, 1] | 0.268 [0.00985, 0.864] | 0 [0, 0] |
| 10028 | power_range09_11 | 6.76e-07 [2.94e-11, 0.000912] | 1.46e-15 [3.65e-22, 3.21e-10] | 2.09e-275 [1.82e-302, 1.27e-249] | 0 [0, 0] |
| 10028 | power_bispectrum_fixed | 1 [1, 1] | 1 [1, 1] | 0.503 [0.0439, 0.957] | 0 [0, 0] |
| 10028 | power_bispectrum_range09_11 | 4.61e-05 [9.26e-09, 0.0141] | 2.11e-12 [2.69e-18, 9.19e-08] | 4.2e-255 [3.57e-281, 2.6e-230] | 0 [0, 0] |
| 10049 | power_fixed | 1 [1, 1] | 1 [1, 1] | 6.82e-21 [1.55e-28, 1.63e-14] | 0 [0, 0] |
| 10049 | power_range09_11 | 1 [1, 1] | 0.994 [0.79, 1] | 1.14e-107 [9.02e-125, 7.56e-92] | 0 [0, 0] |
| 10049 | power_bispectrum_fixed | 1 [1, 1] | 1 [1, 1] | 6.5e-15 [2.23e-21, 1.05e-09] | 0 [0, 0] |
| 10049 | power_bispectrum_range09_11 | 1 [1, 1] | 1 [0.98, 1] | 4.48e-103 [8.42e-120, 1.26e-87] | 0 [0, 0] |

The full CSV retains all 192 projections at n=1,000, 10,000 and 100,000, including both correct-map and removed-map null controls. Maximum point projections of null rejection are 0.0428997 for the true-map control and 0.0434799 for removed-map amplitude-one data; their complete intervals remain in the CSV. These are projections conditional on the realized calibration, not empirical unconditional false-positive rates.

All eight independent first-view physical-cell sums agree with the NUFFT operator within 6.27e-12 per transferred Fourier coordinate. Every saved event count, critical value and first-view score is independently replayed in the accompanying verification record. A threshold perturbation of plus/minus 1e-8 changes 1 count across 64 saved score arrays; the smallest score-to-threshold distance is 8.6e-09. Numerical agreement is not outward-rounded arithmetic.

## Interpretation

The viewing-density bound is an added modeling assumption, not a solution to the failed arbitrary-view certificate. Increasing kappa changes the critical value while the held-out alternative remains Haar; no non-Haar test sample is generated. The resulting power loss measures conservatism of this probability bound, not intrinsic sensitivity to preferred views. Power-only comparisons and all amplitude-range failures remain visible. A score fitted with knowledge of the oracle alternative does not establish practical sensitivity to unknown structural errors. Only one CTF/noise profile and the full region-removal alternative are tested; smaller deviations, unknown noise, image windows and experimentally calibrated viewing/amplitude laws remain unresolved. The construction combines classical binomial inference and nuisance-event domination. Its empirical usefulness and novelty must be assessed separately from its conditional mathematical validity. All three full Fable reviews remain rejections; no new full acceptance review occurred.
