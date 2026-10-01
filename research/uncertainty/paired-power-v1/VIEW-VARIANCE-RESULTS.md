# Conditional-noise and viewing-variance refinement: complete results

This is a selected post-outcome follow-up of all four nonzero ranged-amplitude directions on 10028/10049. It is not a new three-stack positive result; the fixed-amplitude and zero 10076 directions remain in the previous inventory. Scores and thresholds are unchanged. Each stack has 32,768 new Haar calibration views, two groups of 32 independent noises per view, and 131,072 new independent Haar held-out views. Both true and removed candidates and all three matched-simulation-budget procedures are retained.

The three methods use the individual-noise amplitude envelope, grouped amplitude envelope, and grouped viewing-variance bound, respectively. All use the same assumed [.9,1.1] amplitude interval, known CTF/proper Gaussian noise, and supplied viewing-density ratio. The error argument uses independent views as its sample units. Conditional replicas are not counted as independent views. No experimental parameter is calibrated here.

## Removed-candidate rejection projections at n=10,000

Each entry is a binomial projection on the Haar true-map alternative with a pointwise 95% interval from held-out event frequencies. These are not observed repeated-dataset power or simultaneous confidence intervals. Changing kappa changes the null guarantee and threshold, not the held-out viewing distribution.

| Stack | Score | Procedure | kappa=1 | kappa=1.1 | kappa=2 | kappa=5 |
|---|---|---|---:|---:|---:|---:|
| 10028 | power_range09_11 | individual_ratio | 0.0224 [0.00541, 0.0716] | 1.85e-31 [2.81e-34, 9.08e-29] | 0 [0, 0] | 0 [0, 0] |
| 10028 | power_range09_11 | grouped_ratio | 0.0564 [0.0167, 0.148] | 2.9e-29 [5.58e-32, 1.12e-26] | 0 [0, 0] | 0 [0, 0] |
| 10028 | power_range09_11 | view_variance | 0.0499 [0.0143, 0.135] | 0.00108 [0.000154, 0.0058] | 3.48e-10 [9.77e-12, 9.3e-09] | 4.64e-27 [1.15e-29, 1.41e-24] |
| 10028 | power_bispectrum_range09_11 | individual_ratio | 0.0748 [0.0237, 0.184] | 5.52e-28 [1.22e-30, 1.86e-25] | 0 [0, 0] | 0 [0, 0] |
| 10028 | power_bispectrum_range09_11 | grouped_ratio | 0.195 [0.0804, 0.375] | 5.14e-25 [1.61e-27, 1.23e-22] | 0 [0, 0] | 0 [0, 0] |
| 10028 | power_bispectrum_range09_11 | view_variance | 0.184 [0.0746, 0.36] | 0.00816 [0.00162, 0.0315] | 8.38e-09 [3.1e-10, 1.7e-07] | 1.47e-25 [4.32e-28, 3.75e-23] |
| 10049 | power_range09_11 | individual_ratio | 0.505 [0.299, 0.711] | 5.86e-18 [4.54e-20, 5.65e-16] | 0 [0, 0] | 0 [0, 0] |
| 10049 | power_range09_11 | grouped_ratio | 0.561 [0.349, 0.757] | 2.33e-17 [1.96e-19, 2.06e-15] | 0 [0, 0] | 0 [0, 0] |
| 10049 | power_range09_11 | view_variance | 0.546 [0.334, 0.744] | 0.116 [0.041, 0.256] | 3.22e-05 [2.83e-06, 0.000278] | 3.4e-16 [3.41e-18, 2.53e-14] |
| 10049 | power_bispectrum_range09_11 | individual_ratio | 0.591 [0.377, 0.78] | 2.97e-18 [2.22e-20, 2.97e-16] | 0 [0, 0] | 0 [0, 0] |
| 10049 | power_bispectrum_range09_11 | grouped_ratio | 0.702 [0.496, 0.858] | 4.65e-17 [4.12e-19, 3.92e-15] | 0 [0, 0] | 0 [0, 0] |
| 10049 | power_bispectrum_range09_11 | view_variance | 0.681 [0.472, 0.844] | 0.185 [0.0752, 0.362] | 5.1e-05 [4.75e-06, 0.000415] | 7.7e-17 [7.05e-19, 6.29e-15] |

The full CSV retains all 360 projections at n=1,000/10,000/100,000 and kappa=1/1.01/1.1/2/5, including both correct-null controls. The two runs take 391.94 and 408.08 seconds in overlapping Mac processes; their sum is not elapsed wall time. Eight first-view physical-cell sums agree within 8.9e-12 per coordinate. Independent replay checks 120 probability bounds, 360 critical values, all eight held-out count arrays and 512 first-view noise polynomials across 16 amplitude cells.

## Interpretation and limitations

The grouping reduces noise-dependent amplitude selection in the simulation envelope. The paired-group product estimates variation of its conditional mean across views, enabling a less crude density-ratio penalty. Its ingredients are classical empirical Bernstein concentration, conditional covariance and chi-squared/Cauchy--Schwarz robustness, with primary reading scopes documented separately. A numerical improvement is not a novelty verdict. The score still knows an oracle alternative, only one transfer profile and full-removal alternative are tested, and all viewing laws here are simulated Haar. Experimental noise/view/amplitude calibration, heterogeneity, smaller changes and useful learned candidates remain unresolved. All three full Fable reviews remain rejections.
