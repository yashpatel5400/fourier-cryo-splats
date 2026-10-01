# Equal-noise-budget allocation: complete results

Each allocation uses 2,097,152 conditional noise draws per stack: M=32,768 views with two groups of 32 versus M=8,192 with two groups of 128. The larger-group calibration is new; the original one is frozen. Both are evaluated on the same 131,072 fresh test views, distinct from the previously viewed Fisher study. Directions and thresholds are unchanged. This is equal noise-draw count, not equal Fourier calls or wall time, and is not an optimized allocation or repeated-calibration study.

## Mean-envelope and outer-uncertainty tradeoff

| Stack | Score | L | M | Grouped mean | Joint mean upper | Viewing variance upper |
|---|---|---:|---:|---:|---:|---:|
| 10028 | power_matched | 32 | 32768 | 0.516602 | 0.521196 | 0.026187 |
| 10028 | power_matched | 128 | 8192 | 0.513680 | 0.523836 | 0.028093 |
| 10028 | power_fisher | 32 | 32768 | 0.514603 | 0.519054 | 0.024050 |
| 10028 | power_fisher | 128 | 8192 | 0.511803 | 0.521655 | 0.025925 |
| 10028 | power_bispectrum_matched | 32 | 32768 | 0.542132 | 0.546575 | 0.024220 |
| 10028 | power_bispectrum_matched | 128 | 8192 | 0.538295 | 0.548133 | 0.026342 |
| 10028 | power_bispectrum_fisher | 32 | 32768 | 0.522783 | 0.527256 | 0.024466 |
| 10028 | power_bispectrum_fisher | 128 | 8192 | 0.519328 | 0.529238 | 0.026440 |
| 10049 | power_matched | 32 | 32768 | 0.521533 | 0.524328 | 0.005361 |
| 10049 | power_matched | 128 | 8192 | 0.518525 | 0.524590 | 0.006396 |
| 10049 | power_fisher | 32 | 32768 | 0.520390 | 0.523156 | 0.005140 |
| 10049 | power_fisher | 128 | 8192 | 0.517363 | 0.523337 | 0.006075 |
| 10049 | power_bispectrum_matched | 32 | 32768 | 0.523845 | 0.526639 | 0.005350 |
| 10049 | power_bispectrum_matched | 128 | 8192 | 0.520387 | 0.526462 | 0.006432 |
| 10049 | power_bispectrum_fisher | 32 | 32768 | 0.521477 | 0.524240 | 0.005097 |
| 10049 | power_bispectrum_fisher | 128 | 8192 | 0.518385 | 0.524378 | 0.006143 |
| 10076 | power_matched | 32 | 32768 | 0.522836 | 0.525065 | 0.001312 |
| 10076 | power_matched | 128 | 8192 | 0.519811 | 0.524301 | 0.002216 |
| 10076 | power_fisher | 32 | 32768 | 0.522030 | 0.524260 | 0.001325 |
| 10076 | power_fisher | 128 | 8192 | 0.519173 | 0.523671 | 0.002216 |
| 10076 | power_bispectrum_matched | 32 | 32768 | 0.524601 | 0.526833 | 0.001332 |
| 10076 | power_bispectrum_matched | 128 | 8192 | 0.521399 | 0.525901 | 0.002265 |
| 10076 | power_bispectrum_fisher | 32 | 32768 | 0.523249 | 0.525481 | 0.001327 |
| 10076 | power_bispectrum_fisher | 128 | 8192 | 0.520059 | 0.524564 | 0.002262 |

## Fixed-score power comparison

Entries give L=32 / L=128 rejection projections at n=10,000 and amplitude one. The full CSV retains all seven methods, sample sizes, caps, amplitudes and correct-null controls. No cross-allocation minimum is used as a valid test.

| Stack | Score | Method | kappa | 25% deletion | 50% | Full removal |
|---|---|---|---:|---:|---:|---:|
| 10028 | power_matched | view_variance | 1 | 5.294e-05 / 5.48e-06 | 0.1744 / 0.07247 | 1 / 1 |
| 10028 | power_matched | view_variance | 1.1 | 1.166e-45 / 1.645e-50 | 1.997e-29 / 2.605e-33 | 1.137e-06 / 1.942e-08 |
| 10028 | power_matched | cvar_split | 1 | 3.492e-05 / 2.116e-07 | 0.1499 / 0.01705 | 1 / 1 |
| 10028 | power_matched | cvar_split | 1.1 | 2.499e-21 / 2.444e-24 | 4.815e-11 / 3.591e-13 | 0.5052 / 0.2445 |
| 10028 | power_fisher | view_variance | 1 | 7.209e-05 / 7.753e-06 | 0.2295 / 0.1036 | 1 / 1 |
| 10028 | power_fisher | view_variance | 1.1 | 1.285e-42 / 3.455e-48 | 1.866e-26 / 7.769e-31 | 2.412e-05 / 3.291e-07 |
| 10028 | power_fisher | cvar_split | 1 | 4.409e-05 / 2.071e-07 | 0.1947 / 0.02261 | 1 / 1 |
| 10028 | power_fisher | cvar_split | 1.1 | 1.607e-20 / 1.479e-23 | 3.769e-10 / 3.047e-12 | 0.6446 / 0.3617 |
| 10028 | power_bispectrum_matched | view_variance | 1 | 5.15e-06 / 1.11e-06 | 0.04465 / 0.02165 | 0.9999 / 0.9998 |
| 10028 | power_bispectrum_matched | view_variance | 1.1 | 8.608e-47 / 1.33e-51 | 1.714e-31 / 1.965e-35 | 4.206e-10 / 2.774e-12 |
| 10028 | power_bispectrum_matched | cvar_split | 1 | 4.274e-06 / 2.165e-07 | 0.04099 / 0.009572 | 0.9999 / 0.9992 |
| 10028 | power_bispectrum_matched | cvar_split | 1.1 | 1.054e-23 / 1.354e-24 | 1.913e-13 / 4.207e-14 | 0.03833 / 0.02425 |
| 10028 | power_bispectrum_fisher | view_variance | 1 | 4.506e-05 / 7.927e-06 | 0.318 / 0.191 | 1 / 1 |
| 10028 | power_bispectrum_fisher | view_variance | 1.1 | 7.186e-44 / 9.306e-49 | 1.134e-25 / 2.126e-29 | 0.0002901 / 1.163e-05 |
| 10028 | power_bispectrum_fisher | cvar_split | 1 | 3.222e-05 / 5.369e-07 | 0.29 / 0.07552 | 1 / 1 |
| 10028 | power_bispectrum_fisher | cvar_split | 1.1 | 1.995e-21 / 1.787e-23 | 1.017e-09 / 4.647e-11 | 0.8378 / 0.692 |
| 10049 | power_matched | view_variance | 1 | 0.0176 / 0.01515 | 0.8387 / 0.8236 | 1 / 1 |
| 10049 | power_matched | view_variance | 1.1 | 7.969e-12 / 2.578e-13 | 0.0001325 / 1.816e-05 | 0.9332 / 0.8454 |
| 10049 | power_matched | cvar_split | 1 | 0.01054 / 0.003826 | 0.7849 / 0.6656 | 1 / 1 |
| 10049 | power_matched | cvar_split | 1.1 | 3.443e-08 / 2.681e-07 | 0.01069 / 0.02745 | 0.9978 / 0.9994 |
| 10049 | power_fisher | view_variance | 1 | 0.01327 / 0.01196 | 0.8041 / 0.7928 | 1 / 1 |
| 10049 | power_fisher | view_variance | 1.1 | 6.433e-12 / 3.197e-13 | 0.0001084 / 1.888e-05 | 0.9544 / 0.8972 |
| 10049 | power_fisher | cvar_split | 1 | 0.007789 / 0.003079 | 0.7441 / 0.6313 | 1 / 1 |
| 10049 | power_fisher | cvar_split | 1.1 | 2.581e-08 / 2.803e-07 | 0.00882 / 0.02671 | 0.9987 / 0.9997 |
| 10049 | power_bispectrum_matched | view_variance | 1 | 0.01466 / 0.01622 | 0.8204 / 0.8307 | 1 / 1 |
| 10049 | power_bispectrum_matched | view_variance | 1.1 | 4.758e-12 / 2.687e-13 | 9.907e-05 / 1.865e-05 | 0.9256 / 0.851 |
| 10049 | power_bispectrum_matched | cvar_split | 1 | 0.008209 / 0.004397 | 0.7569 / 0.6829 | 1 / 1 |
| 10049 | power_bispectrum_matched | cvar_split | 1.1 | 2.546e-08 / 4.191e-07 | 0.009271 / 0.03345 | 0.9975 / 0.9996 |
| 10049 | power_bispectrum_fisher | view_variance | 1 | 0.02175 / 0.02073 | 0.8469 / 0.8421 | 1 / 1 |
| 10049 | power_bispectrum_fisher | view_variance | 1.1 | 2.849e-11 / 1.163e-12 | 0.0002216 / 3.515e-05 | 0.9688 / 0.9192 |
| 10049 | power_bispectrum_fisher | cvar_split | 1 | 0.01193 / 0.005872 | 0.783 / 0.6991 | 1 / 1 |
| 10049 | power_bispectrum_fisher | cvar_split | 1.1 | 6.942e-08 / 9.433e-07 | 0.01297 / 0.04224 | 0.9992 / 0.9999 |
| 10076 | power_matched | view_variance | 1 | 0.0004687 / 0.0007662 | 0.005552 / 0.008216 | 0.1963 / 0.2375 |
| 10076 | power_matched | view_variance | 1.1 | 9.916e-09 / 4.261e-10 | 6.345e-07 / 4.024e-08 | 0.000787 / 0.0001155 |
| 10076 | power_matched | cvar_split | 1 | 0.0001933 / 0.0003775 | 0.002719 / 0.004668 | 0.1367 / 0.1801 |
| 10076 | power_matched | cvar_split | 1.1 | 1.305e-09 / 7.489e-07 | 1.077e-07 / 2.644e-05 | 0.0002311 / 0.009179 |
| 10076 | power_fisher | view_variance | 1 | 0.0002797 / 0.0004336 | 0.004801 / 0.006762 | 0.1818 / 0.2153 |
| 10076 | power_fisher | view_variance | 1.1 | 4.33e-09 / 1.337e-10 | 4.923e-07 / 2.432e-08 | 0.0006545 / 7.975e-05 |
| 10076 | power_fisher | cvar_split | 1 | 0.0001032 / 0.000178 | 0.002184 / 0.003364 | 0.1212 / 0.1518 |
| 10076 | power_fisher | cvar_split | 1.1 | 4.798e-10 / 2.674e-07 | 7.375e-08 / 1.642e-05 | 0.000175 / 0.00672 |
| 10076 | power_bispectrum_matched | view_variance | 1 | 0.0004578 / 0.0008595 | 0.007028 / 0.01144 | 0.1909 / 0.2438 |
| 10076 | power_bispectrum_matched | view_variance | 1.1 | 9.518e-09 / 4.071e-10 | 9.607e-07 / 6.348e-08 | 0.0007333 / 0.0001065 |
| 10076 | power_bispectrum_matched | cvar_split | 1 | 0.0002034 / 0.0004578 | 0.003723 / 0.007028 | 0.1367 / 0.1909 |
| 10076 | power_bispectrum_matched | cvar_split | 1.1 | 1.412e-09 / 9.763e-07 | 1.862e-07 / 4.828e-05 | 0.0002306 / 0.01021 |
| 10076 | power_bispectrum_fisher | view_variance | 1 | 0.0003731 / 0.000707 | 0.005492 / 0.009066 | 0.1769 / 0.2276 |
| 10076 | power_bispectrum_fisher | view_variance | 1.1 | 6.858e-09 / 2.853e-10 | 6.221e-07 / 3.937e-08 | 0.0006118 / 8.665e-05 |
| 10076 | power_bispectrum_fisher | cvar_split | 1 | 0.0001517 / 0.0003223 | 0.002687 / 0.004893 | 0.1214 / 0.1667 |
| 10076 | power_bispectrum_fisher | cvar_split | 1.1 | 8.838e-10 / 6.023e-07 | 1.055e-07 / 2.833e-05 | 0.0001755 / 0.007958 |

## Controls and verification

Largest correct-null n=10,000 projection is 7.8957418e-05 [7.7946547e-06,0.00060523612] at 10076, power_bispectrum_fisher, L=128, grouped_ratio, kappa=1.0, amplitude=1.1.
Largest correct-null group count is 1/128 [0.00019778,0.042759]. Both calibration realizations are fixed. All cells share draws and cannot be pooled as independent replications.

Independent replay checks 12 unchanged score designs, 420 new bounds, 3360 minimal critical values, 37800 projections and 12600 group vectors. It checks 3072 first-view cubics plus direct physical scores. It does not replay every physical view or provide a certified rounding bound.

The declared known-simulator assumptions and whole-candidate interpretation remain unchanged. This comparison does not establish experimental noise/view calibration, local occupancy coverage, an optimal allocation or a new full-review verdict.
