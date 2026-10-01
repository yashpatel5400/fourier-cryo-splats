# Saved-array alignment diagnosis after review 3

Post hoc analysis of all 600 frozen trials. All numerical inputs are hash-checked; no new random observations are generated. The true-pose controls are newly fitted deterministic operators on the saved images. The reduced-radius rows retain their original weights and explicitly flag whether their class contains the generator. All individual Monte Carlo intervals are nonsimultaneous.

## Audit-weight results

| Stack | Alignment | Target | Mean error / SD, same / independent | Noise-only coverage, same / independent | RMSE / pilot error, same / independent |
|---|---|---|---|---|---|
| 10028 | true_pose | center | -0.367 / -0.400 | 0.920 / 0.925 | 0.061 / 0.062 |
| 10028 | true_pose | contrast | -0.099 / -0.139 | 0.965 / 0.965 | 0.276 / 0.311 |
| 10028 | oracle_reference | center | -0.258 / -0.539 | 0.925 / 0.895 | 0.056 / 0.069 |
| 10028 | oracle_reference | contrast | -0.156 / 0.137 | 0.960 / 0.935 | 0.265 / 0.297 |
| 10028 | independent_pilot | center | -1.177 / -1.266 | 0.815 / 0.750 | 0.085 / 0.092 |
| 10028 | independent_pilot | contrast | 0.217 / 0.239 | 0.965 / 0.940 | 0.299 / 0.325 |
| 10049 | true_pose | center | 0.081 / 0.070 | 0.930 / 0.980 | 1.432 / 1.391 |
| 10049 | true_pose | contrast | 0.198 / 0.194 | 0.930 / 0.955 | 1.029 / 0.999 |
| 10049 | oracle_reference | center | 0.737 / -0.295 | 0.900 / 0.925 | 1.590 / 1.478 |
| 10049 | oracle_reference | contrast | -0.128 / 0.656 | 0.950 / 0.895 | 0.757 / 1.018 |
| 10049 | independent_pilot | center | 0.638 / -0.557 | 0.905 / 0.925 | 1.600 / 1.594 |
| 10049 | independent_pilot | contrast | 1.007 / 1.349 | 0.860 / 0.720 | 1.127 / 1.350 |
| 10076 | true_pose | center | -0.282 / -0.173 | 0.920 / 0.960 | 0.343 / 0.314 |
| 10076 | true_pose | contrast | -0.068 / -0.161 | 0.955 / 0.925 | 0.758 / 0.872 |
| 10076 | oracle_reference | center | 0.363 / -0.384 | 0.955 / 0.915 | 0.317 / 0.330 |
| 10076 | oracle_reference | contrast | 0.165 / -0.024 | 0.950 / 0.945 | 0.689 / 0.797 |
| 10076 | independent_pilot | center | 0.941 / -0.325 | 0.860 / 0.935 | 0.433 / 0.332 |
| 10076 | independent_pilot | contrast | -0.162 / 0.081 | 0.965 / 0.950 | 0.718 / 0.792 |

## Continuous class at realized poses

These are oracle diagnostics at the realized pose pair, not deployable intervals or uniform pose-set bounds. All norms use the continuous cube; Gauss quadrature errors are bounded analytically, while NUFFT and rounding are only checked numerically.

| Stack | Template | Target | Total bias / no data, median | Pose-only bias / no data, median | First-order half-width / no data, median |
|---|---|---|---|---|---|
| 10028 | oracle_reference | center | 0.364 | 0.355 | 4.527 |
| 10028 | oracle_reference | contrast | 0.415 | 0.404 | 5.411 |
| 10028 | independent_pilot | center | 0.469 | 0.462 | 4.857 |
| 10028 | independent_pilot | contrast | 0.526 | 0.517 | 5.887 |
| 10049 | oracle_reference | center | 0.464 | 0.456 | 5.052 |
| 10049 | oracle_reference | contrast | 0.487 | 0.474 | 5.307 |
| 10049 | independent_pilot | center | 0.503 | 0.495 | 5.238 |
| 10049 | independent_pilot | contrast | 0.532 | 0.519 | 5.472 |
| 10076 | oracle_reference | center | 0.306 | 0.295 | 3.215 |
| 10076 | oracle_reference | contrast | 0.341 | 0.325 | 3.849 |
| 10076 | independent_pilot | center | 0.326 | 0.315 | 3.365 |
| 10076 | independent_pilot | contrast | 0.374 | 0.359 | 4.024 |

## Exact pilot-to-generator distance

- 10028: 1.262007; all 6 true-pose fits converged: True; maximum replay discrepancy 0.
- 10049: 1.045950; all 6 true-pose fits converged: True; maximum replay discrepancy 0.
- 10076: 1.027680; all 6 true-pose fits converged: True; maximum replay discrepancy 0.

The machine-readable summaries retain all 156 estimator/image cells, 78 paired comparisons and 2,400 realized envelope calculations. Class-radius and noise-only coverage are diagnostic, not evidence that these smaller classes are valid experimentally.
