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
