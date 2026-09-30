# Completed centered-noise sensitivity

30 September 2026. All twelve unchanged fixed estimators, all 48 pose-class
intervals, and all three unchanged cubic estimators completed the declared
centered-noise procedure. All uncentered SD bounds replayed against their
published values. The contrast uses 127 independent rows from 128 calibration
exposures under the common Gaussian covariance assumption; unequal signal
means are permitted. These pixels had already been analyzed.

For the twelve fixed estimators, centered/uncentered SD ratios range from
0.43739 to 0.82343. Zero exclusions change from six to ten at fixed poses and
from four to six for shift-only sensitivity. They remain zero for one- and
two-degree rotations. The four one-degree and twelve two-degree no-data
fallbacks remain. Five registered reference targets now lie outside: 10076
regions 1 and 3 at fixed/shift-only poses, and its matched center at fixed poses.
Thus 43 of 48 fixed-weight intervals contain the registered reference; the
three cubic intervals also contain it. These are agreements/disagreements,
not coverage labels: 10076 is heterogeneous and EMD-8434 is one Class A state. Centering does not establish
pose radii, common covariance, homogeneity, or a physical density class.

| Completed cubic estimator | Center | Centered SD upper | Half-width | Fraction of no data |
| --- | ---: | ---: | ---: | ---: |
| Original | 2.46178 | .55294 | 3.07893 | .25323 |
| Coordinate metric | 2.48844 | .55466 | 3.60408 | .29642 |
| Reduced subspace | 2.46178 | .55294 | 3.07784 | .25314 |

None excludes zero or falls back to no data; all contain the registered
reference target 2.48846. The original estimator's SD upper falls from 1.01002
to .55294, still above its supplied simulation SD .22191. This indicates that
the uncentered bound included a large removable common component; it does not
identify the remaining contribution as pure noise. The exposure-dependence
inflation is unchanged. No minimum of the two noise bounds is selected.

The three new contrast tests pass, including a separately generated unequal-
mean, correlated Gaussian-covariance check. The full suite at this checkpoint
reports 230 passed, 1 skipped and 1 intentional warning in 90.75 seconds.
This precedes the later cubic enrichment module, whose two separate prerequisite
tests pass; no combined full-suite count is inferred. Tests are implementation
checks, not experimental coverage evidence.
