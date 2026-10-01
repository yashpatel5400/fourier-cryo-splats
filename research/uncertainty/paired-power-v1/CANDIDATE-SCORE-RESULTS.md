# Candidate-only local-perturbation study: complete results

All three previously fitted Gaussian maps are used, with no external reference map for region or score design. Each candidate supplies a 20-Angstrom smoothed peak and its own prescribed 25% local-deletion direction. Four scores compare power/combined moments and raw/scale-orthogonal directions. Training uses 4,096 fresh Haar views, calibration 32,768 views with two groups of 32 noises, and testing 131,072 independent Haar views per stack. All noise, transfer and viewing assumptions remain specified simulators. These are new simulations of fitted candidates, not new reconstructions or applications to experimental test particles.

## Paired-variance projections at n=10,000 and amplitude one

Each entry is a binomial projection, not observed repeated-dataset power. The complete CSV includes pointwise 95% intervals, both other amplitudes, every sample size and all seven calibration procedures. Rounded tiny values are displayed in scientific notation. Exact non-rejection and numerical underflow are distinguished by the saved critical counts.

| Stack | Score | kappa | 10% removal | 25% | 50% | 100% |
|---|---|---:|---:|---:|---:|---:|
| 10028 | power_raw | 1 | <1e-300 | <1e-300 | 1.7e-261 | 1.73e-143 |
| 10028 | power_raw | 1.1 | <1e-300 | <1e-300 | <1e-300 | 3.83e-238 |
| 10028 | power_scale_orthogonal | 1 | 1.41e-09 | 2.25e-05 | 0.178 | 1 |
| 10028 | power_scale_orthogonal | 1.1 | 1.09e-59 | 2.82e-47 | 1.27e-29 | 4.68e-07 |
| 10028 | power_bispectrum_raw | 1 | 7.66e-268 | 2.19e-233 | 2.59e-182 | 6.9e-100 |
| 10028 | power_bispectrum_raw | 1.1 | <1e-300 | <1e-300 | 2.07e-305 | 5.61e-196 |
| 10028 | power_bispectrum_scale_orthogonal | 1 | 1.19e-10 | 8.93e-07 | 0.0169 | 1 |
| 10028 | power_bispectrum_scale_orthogonal | 1.1 | 1.49e-59 | 4.37e-49 | 1.15e-33 | 3.95e-11 |
| 10049 | power_raw | 1 | 2.08e-45 | 1.81e-24 | 3.21e-05 | 1 |
| 10049 | power_raw | 1.1 | 3.42e-93 | 1.72e-61 | 1.51e-25 | 0.759 |
| 10049 | power_scale_orthogonal | 1 | 0.000403 | 0.0776 | 0.926 | 1 |
| 10049 | power_scale_orthogonal | 1.1 | 4.39e-16 | 4.82e-10 | 0.000574 | 0.983 |
| 10049 | power_bispectrum_raw | 1 | 1.17e-43 | 1.51e-22 | 0.000483 | 1 |
| 10049 | power_bispectrum_raw | 1.1 | 8.87e-91 | 1.8e-58 | 1.44e-22 | 0.971 |
| 10049 | power_bispectrum_scale_orthogonal | 1 | 0.00019 | 0.0625 | 0.936 | 1 |
| 10049 | power_bispectrum_scale_orthogonal | 1.1 | 6.9e-17 | 2.05e-10 | 0.000696 | 0.984 |
| 10076 | power_raw | 1 | 1.65e-17 | 1.53e-14 | 1.34e-10 | 3.51e-05 |
| 10076 | power_raw | 1.1 | 1.96e-26 | 1.06e-22 | 1.38e-17 | 4.72e-10 |
| 10076 | power_scale_orthogonal | 1 | 0.000458 | 0.00263 | 0.0224 | 0.342 |
| 10076 | power_scale_orthogonal | 1.1 | 7.55e-09 | 1.41e-07 | 6.77e-06 | 0.00296 |
| 10076 | power_bispectrum_raw | 1 | 5.04e-19 | 1.03e-15 | 1.55e-11 | 1.6e-05 |
| 10076 | power_bispectrum_raw | 1.1 | 9.58e-29 | 1.4e-24 | 3.6e-19 | 8.22e-11 |
| 10076 | power_bispectrum_scale_orthogonal | 1 | 0.000107 | 0.000846 | 0.013 | 0.246 |
| 10076 | power_bispectrum_scale_orthogonal | 1.1 | 5.71e-10 | 1.64e-08 | 1.99e-06 | 0.00106 |

## Classical comparison on full deletion, kappa=1.1, n=10,000

These separately declared tests share draws. No minimum across methods is used as a valid combined test. The table is restricted by score type (the scale-orthogonal arm); all raw-arm results remain in the CSV.

| Stack | Score | Paired variance | CVaR split | CVaR DKW | Grouped ratio | Unpaired variance |
|---|---|---:|---:|---:|---:|---:|
| 10028 | power_scale_orthogonal | 4.68e-07 | 0.426 | 0.22 | 5.19e-07 | 3.09e-07 |
| 10028 | power_bispectrum_scale_orthogonal | 3.95e-11 | 0.0219 | 0.00267 | 3.18e-14 | 3.63e-13 |
| 10049 | power_scale_orthogonal | 0.983 | 1 | 0.998 | 0.000131 | 0.781 |
| 10049 | power_bispectrum_scale_orthogonal | 0.984 | 1 | 0.998 | 0.000118 | 0.786 |
| 10076 | power_scale_orthogonal | 0.00296 | 0.00113 | 9.35e-05 | 5.6e-28 | 4.77e-07 |
| 10076 | power_bispectrum_scale_orthogonal | 0.00106 | 0.000459 | 2.7e-05 | 1.52e-29 | 9.83e-08 |

## Correct-null controls, scope and verification

Across every method/score/kappa/amplitude at n=10,000, the largest correct-null point projection is 0.0537238. Across all 1260 correct-null repeated-group cells, the largest observed count is 7/128 (pointwise exact 95% interval [0.022267, 0.10943]). Calibration is fixed, cells share draws, and these cells must not be pooled as independent repetitions. This is not a calibration-repeat or unconditional-error experiment.

The full report retains 18,900 scalar projections (1,260 method rows) and 6,300 actual group-outcome cells. Both smaller deletions and all zero-deletion/amplitude controls remain visible. The removed-region L2 fractions are 10028: 0.08711, 10049: 0.26546, 10076: 0.14639; multiplying these by the deletion amount gives the relative L2 change in the candidate. Normalization and the transfer/noise fixture are specified, not estimated experimental signal-to-noise ratios.

Independent FFT convolution reproduces all three candidate-only region centers. Separate feature arithmetic and least-squares projection reproduce all 12 directions within 7.87e-15. The check includes 768 first-view calibration cubics, 180 direct held-out scores, every held-out event count, 1680 minimal critical values, all 6,300 group vectors and all 18,900 rejection projections. Maximum direct-score discrepancy is 2.45e-11. It is a numerical replay, not a global floating-point certificate.

The three Mac processes take 10028: 408.63s, 10049: 404.70s, 10076: 403.02s. They overlap; summed process durations are not elapsed wall time.

## Scientific interpretation

Candidate-driven design removes the need to know an external true-map alternative for training this score. It does not remove the hand-specified deletion family or calibrate experimental imaging nuisance. Scale orthogonalization improves several full/half-deletion cases, while 10--25% changes remain weak. Stronger view-law protection again removes much of the sensitivity. Unlike the previous selected oracle contrasts, these candidate scores show substantial empirical between-view variation; the paired variance decomposition is retained for every score. On 10028 and 10049 the CVaR comparator can be tighter than the paired variance method. This reversal is retained and rules out a general empirical superiority claim from the earlier experiment.

A rejection concerns the complete candidate under the stated simulator. Other structural/imaging errors can change the same moment score, so neither rejection nor non-rejection is a unique local-error attribution or confidence interval for regional occupancy. All three full Fable reviews still reject. No new full acceptance assessment has occurred.
