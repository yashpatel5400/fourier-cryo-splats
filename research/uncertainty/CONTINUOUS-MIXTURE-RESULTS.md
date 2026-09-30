# Continuous-orientation computation: unresolved global bounds

30 September 2026. All three prospectively specified Gaussian-pilot simulations
completed. These are computational diagnostics on 128 independent synthetic
images per acquisition geometry, with supplied Gaussian noise, fixed CTFs and
zero shifts. They do not measure experimental coverage or structural-test power.
The truth orientations are used only in the separately reported feasible oracle
support and local-bound diagnostic, never as input to the global solver.

Each initial global computation exhausted 8,192 splits, leaving 8,224 boxes
covering SO(3) and 16,416 evaluated centers. None reached the one-log-unit gap
target. The three solves took 430.75, 379.51 and 390.67 seconds; total runner
time was 1,208.25 seconds including data construction and saving.

| Geometry | Center-mixture lower | Initial global upper | Upper after envelope-weight fitting | Gap above best recorded feasible value |
|---|---:|---:|---:|---:|
| 10028 | -82017.232 | -55185.096 | -60233.542 | 20253.222 |
| 10049 | -81391.289 | -60843.210 | -62998.213 | 17988.167 |
| 10076 | -80537.454 | -60663.177 | -63460.256 | 17077.197 |

All likelihoods include the same Gaussian normalization, -51754.618.
The best recorded feasible value additionally considers a uniform mixture on
the 128 generating rotations, evaluated after fitting. That mixture is a
feasible viewing measure, not the Haar likelihood or the global maximizer.
Its values are -80486.764, -80986.380 and -80627.370, respectively. Consequently
these are computable upper-to-feasible gaps, not known errors of the upper bound.

The declared anchor follow-up fits weights on the cell envelopes. Its finite
relaxation gaps are 0.029497, 0.009910 and 0.009773. Geometry 10028 reached its
500-iteration limit; the others met the 0.01 tolerance. The best and final
covers happen to coincide. Optimizing envelope weights improves the original
upper values substantially but leaves the large spatial relaxation unresolved.
An envelope-mixture primal is not a feasible continuous likelihood lower bound.
The follow-up's 708.23-second runner time includes waiting for its parent;
per-fit times are retained separately in every case record.

The curvature follow-up completed all 54 prescribed local cases in 8.85 seconds.
Each uses three fixed particle positions per geometry, a box centered at its
true simulated orientation, and nine local optimization starts. The searches
give feasible values, not guaranteed box maxima. Median upper-to-feasible gaps:

| Sum of Euler half-widths | 10028 original / refined | 10049 original / refined | 10076 original / refined |
|---|---:|---:|---:|
| 0.1 degrees | 0.08905 / 0.001658 | 0.08970 / 0.002054 | 0.04644 / 0.000751 |
| 0.5 degrees | 2.336 / 0.2351 | 2.367 / 0.2923 | 1.234 / 0.1058 |
| 1 degree | 9.865 / 2.194 | 10.05 / 2.747 | 5.308 / 0.9827 |
| 2 degrees | 23.44 / 23.26 | 35.04 / 29.84 | 21.59 / 10.61 |
| 5 degrees | 52.12 / 52.12 | 68.31 / 68.31 | 47.74 / 47.74 |
| 10 degrees | 85.71 / 85.71 | 100.5 / 100.5 | 77.34 / 77.34 |

The second-order log-kernel enclosure helps small local boxes but does not
resolve the coarse global cover. These outcomes motivate an independently
declared refinement change; they do not establish a useful continuous test.
All three outcome families retain source/input hashes, stopping statuses,
arrays for replay, and every prescribed case. The arithmetic bounds remain
real-arithmetic results evaluated in ordinary floating point.


## Envelope-guided refinement and reviewer-motivated disks

The separately declared refinement completes on all three geometries. The
first two stop at the 900-second wall budget; the third exhausts 32,768
additional splits. Final saved-anchor replay is included in the reported upper.

| Geometry | Added splits | Leaves | Feasible lower | Upper | Gap | Solver seconds |
|---|---:|---:|---:|---:|---:|---:|
| 10028 | 30528 | 38752 | -80927.973 | -65036.414 | 15891.560 | 900.70 |
| 10049 | 32000 | 40224 | -80990.141 | -66500.405 | 14489.736 | 901.48 |
| 10076 | 32768 | 40992 | -80419.801 | -67521.974 | 12897.827 | 896.49 |

Total runner time is 2,718.27 seconds. None reaches the one-log-unit target.
The final envelope EM fits are also nonconverged, with finite-relaxation gaps
69.12, 52.14 and 89.67; those are distinct from the much larger spatial gaps.
Considering the previously reported oracle-support feasible values changes
the final upper-to-best-feasible gaps to 15,450.350, 14,485.974 and 12,897.827.
The optimizer did not receive those oracle values or rotations. Wall exit in
the first run occurs between scheduled feasible refits, as disclosed in L6 of
the focused review; this can only weaken the retained lower, not the upper.

The focused Fable audit motivated a separately declared per-frequency Taylor
disk diagnostic. All three independent implementation tests pass; all 54
prescribed local boxes and 285 prescribed coarse-cover cells complete in
14.10 seconds, including the fixed-map orientation-contrast evaluation.
No coarse-cell improvement exceeds 1e-9 on any of 36,480 image/cell coordinates;
the three unchanged-anchor cover uppers remain exactly unchanged. At two-degree
local boxes, median gaps change from 23.261 to 22.920 and 29.836 to 27.710 on
10028/10049. The corresponding 10076 median and all other reported budget
medians are unchanged. This local improvement does not solve the global problem.

All three half-plane frequency arrays pass the no-DC/no-duplicate/no-conjugate
check. From 256 independently seeded Haar orientation evaluations per map,
the median best-recorded-center minus 99th-percentile random-orientation log
kernel is 26.995, 22.760 and 7.513, against median envelope-minus-best-center
slacks of 166.454, 141.461 and 128.417. These finite diagnostics support the
coarse-enclosure diagnosis, not an impossibility statement at arbitrary budgets.

Historical anchor labels are clarified without overwriting the files:
`independent_image` is the sum of per-image maximum log envelopes plus Gaussian
normalization; `row_max` is the scaled Lindsay upper *using row maxima as
anchors*. The latter can be smaller, since its scale correction is nonpositive.
The difference between the independent-image upper and the fitted envelope
upper is 289.525, 253.101 and 110.083, small relative to the spatial gaps.
