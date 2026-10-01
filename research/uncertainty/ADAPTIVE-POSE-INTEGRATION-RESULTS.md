# Adaptive orientation integration: first frozen gate

1 October 2026 UTC. All 384 observed-image cases complete: 64 fixed source indices × unchanged/25%-deleted image × three stacks. Every case uses eight optimized proposal modes and two independent 8,192-draw importance banks. The source indices, tolerance and ESS requirements were [declared before outcomes](ADAPTIVE-POSE-INTEGRATION-PROTOCOL.md). All original maps, simulator assumptions and images remain unchanged.

**Overall gate: FAIL.** Passing requires at least 90% of images within .01 log-ratio units between banks and all four median ESS values at least 256 on every stack. No stack is omitted.

| EMPIAR | Images within tolerance | Median ESS: bank0 null / alt / bank1 null / alt | Stack gate | Minutes |
|---|---:|---:|---|---:|
| 10028 | 128/128 (1.0000) | 4523.5 / 4530.4 / 4530.4 / 4521.0 | pass | 17.93 |
| 10049 | 103/128 (0.8047) | 2690.4 / 2675.8 / 2603.9 / 2729.7 | fail | 18.05 |
| 10076 | 79/128 (0.6172) | 226.6 / 191.6 / 246.4 / 243.0 | fail | 18.01 |

Times are concurrent wall times on the M4 Pro, including image-specific physical Fourier evaluations. They cannot be interpreted as speedups over the cached uniform-template matrix calculation, which shares templates across images. No GPU was rented.

## Error and optimizer diagnostics

| EMPIAR | Maximum ratio discrepancy | RMS ratio discrepancy | RMS log-integral discrepancy | RMS last-prefix ratio change | Successful local fits |
|---|---:|---:|---:|---:|---:|
| 10028 | 0.00353397 | 0.000731385 | 0.0150872 | 0.000626213 | 1021/1024 |
| 10049 | 0.0937663 | 0.0199069 | 0.15 | 0.0120924 | 1024/1024 |
| 10076 | 0.0738422 | 0.0183302 | 0.16173 | 0.011786 | 1024/1024 |

The same proposal is shared by an image’s two integration banks. Agreement can miss a common uncovered mode; [the numerical requirements note](NUMERICAL-LIKELIHOOD-REQUIREMENTS.md) gives an elementary counterexample. Effective sample size is a weight-concentration diagnostic, not a coverage guarantee. Neither a successful BFGS status nor a positive local Hessian establishes global pose recovery. All 3,072 local-fit outcomes, including nonconvergence, duplicate modes and any reverted terminal point, are retained. No fitted proposal covariance is presented as a pose confidence region.

## Independent replay

- 10028: all 128 proposal densities and integral summaries replay. Maximum log-density discrepancy 1.6e-14; maximum log-integral discrepancy 5.68e-14; six direct cell-sum residual checks differ by at most 8.25e-11. The prespecified gate verdict reproduces.
- 10049: all 128 proposal densities and integral summaries replay. Maximum log-density discrepancy 2.8e-14; maximum log-integral discrepancy 5.68e-14; six direct cell-sum residual checks differ by at most 3.24e-11. The prespecified gate verdict reproduces.
- 10076: all 128 proposal densities and integral summaries replay. Maximum log-density discrepancy 1.24e-14; maximum log-integral discrepancy 2.84e-14; six direct cell-sum residual checks differ by at most 1.57e-11. The prespecified gate verdict reproduces.

The independent proposal check uses full transformed four-dimensional Gaussian covariance coordinates; it does not call the implemented relative-quaternion density. The physical check explicitly sums all 64³ cells at all retained frequencies for six selected orientations per stack, covering both candidate maps. These checks establish consistency of the computation at those points, not a bound on the remaining orientation integral.

## Decision boundary

The first numerical gate does not by itself create an uncertainty method or address measured nonuniform views, image-specific CTFs, amplitude, colored noise, map error, regional specificity or calibration repetitions. The moment branch remains frozen. Under the review-4 stopping rule, at most one substantive integration revision is permitted, directed at the measured failure, with this same cohort, draw count and gate. A second missed prediction ends this numerical branch. A focused independent methodological consultation will inform whether that revision and a specific statistical estimand are worth pursuing. No new full-paper acceptance verdict is claimed.

All per-image results and importance points are preserved under `results/uncertainty/development/adaptive-pose-integration-v1`; the compact CSV also retains the original uniform-bank scores at the same images.
