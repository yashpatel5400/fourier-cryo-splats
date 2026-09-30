# Original-code CryoLike comparison: both scores and both grids

30 September 2026. All 24 declared cases complete with the unchanged author
source at commit `a413ffd265e2815c9f81c492a3d8bfedaf35d737`. The comparison uses
115, 71 and 128 distinct exposures from the previously examined prediction
cohort, one image per exposure, and no supplied orientation/translation during
scoring. Deposited CTFs and supplied-pose training of the learned maps remain.
This is exploratory scoring, not a new reconstruction or uncertainty guarantee.

The original CPU cross-correlation unit test passed before the study. The local
adapter's asymmetric direct-transform tests verify axis and physical-unit
conventions. The author's Fourier--Bessel scoring and integrated-score kernels
are unchanged. Explicit settings, all selected images, source hashes, template
grids and every per-image score are retained in the case records and arrays.

All cases use the radius-twelve band and both 0.4 and 0.2 radian viewing grids
(173 and 470 templates), 128 in-plane angles, and 25 displacements. Total wall
time is 1,218.87 seconds under concurrent Mac load. Per-case times including
setup range from 25.09 to 81.02 seconds; peak process RSS is 1,514,831,872 bytes.
These timings are not a CUDA benchmark or a fine-resolution scalability claim.

Both grids tell the same qualitative story for maximum correlation. Gaussian
and voxel maps score above the locked neural map on 10028 and 10049 and below
it on 10076. The 18-contrast-family percentile bootstrap intervals exclude zero
for these twelve contrasts. This reverses some earlier fixed-pose predictive
rankings under a different, normalized, searched-pose score; it does not show
that either study is an accuracy ground truth.

For integrated log score, the neural map scores above Gaussian and voxel maps
on 10028 and 10076, with the corresponding family intervals excluding zero.
On 10049 those intervals include zero for both methods and grids. The deposited
map's integrated score is above neural on 10028/10049 and below it on 10076;
all six family intervals exclude zero. Maximum correlation instead puts the
deposited map below neural on 10049 as well as 10076. Both metrics and all 36
contrasts are retained; there is no selected favorable score or grid.

Each interval uses 20,000 paired exposure bootstrap draws. The intervals are
descriptive approximations on reused data, not exact familywise guarantees.
Scores depend on different orientation/amplitude/noise treatments and template
normalization. The author's integrated score is not substituted for a normalized
raw-pixel predictive density or used as an e-value. These results supply an
external structural-scoring comparison, not matched-power evidence for the
proposed continuous likelihood test, whose global denominator remains loose.
