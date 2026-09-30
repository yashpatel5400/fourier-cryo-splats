# A shared unknown noise scale: bounded follow-up

30 September 2026, after the discrete screen and before this follow-up.
The per-image/cell variance relaxation lost the local detections. Test whether
profiling one shared positive scale reduces that cost. Reuse all existing
observations and the oracle numerator; this is post-outcome algorithm
development. Keep all three stacks, all four maps and every one of the 16
repeats. Do not change the discrete 64-view catalog or introduce continuous
pose claims. All old outcomes remain unchanged.

For residual distances r_ic in d real coordinates, each Gaussian kernel is
increasing in sigma below r_ic/sqrt(d) and decreasing above it. Therefore the
joint mixture likelihood cannot attain a larger maximum outside
`[min_ic r_ic/sqrt(d), max_ic r_ic/sqrt(d)]` than inside it. This data-dependent
bracket covers the maximum over all positive sigma pointwise. If the lower
endpoint is zero, retain an unresolved infinite upper bound in this prototype;
do not insert an arbitrary positive cutoff.

Use precision t=1/sigma² and separate the common normalization from
`F(t)=max_w sum_i log(sum_c w_c exp(-r_ic² t/2))`. For every fixed w,
log-sum-exp of affine functions is convex in t. The supremum F is therefore
convex. On an interval [a,b], a chord through upper bounds on F(a) and F(b)
upper-bounds F throughout. Adding `(nd/2) log(t) - (nd/2) log(2pi)` gives a
concave scalar upper curve whose maximum is analytic: if its chord slope s
is negative, clip `-nd/(2s)` to [a,b]; otherwise use b. Obtain the endpoint
bounds with the existing mixture dual. This retains the **common** variance
and avoids independently maximizing the variance in each image's kernel.

Evaluate a feasible mixture at each upper curve's maximizer and at unit scale
clipped to the original bracket. Cache every precision evaluation. Split the
scale interval with largest upper bound in half on the log scale. The maximum
leaf upper bound covers all positive scales. Stop at a one-log-unit global
gap or 65 evaluated intervals. Mixture solves use at most 2,000 iterations
and a 1e-4 gap target. Stopping early retains the upper bound; it cannot license
an optimistic score. Ordinary floating point is not certified interval
arithmetic. Convex chords and branch and bound are established tools; this is
not a new general testing principle.

An initial implementation used separate per-image variance-interval
envelopes. Its one-component numerical test failed the requested stopping
tolerance within 129 nodes, without an observed bound violation. That source,
test and failure log are archived. The precision-chord replacement and this
final protocol precede any empirical shared-scale outcome.

Before data evaluation, check the one-component case against its analytic
variance MLE and a two-component case against an independent dense variance
grid with scalar weight optimization. Verify that terminal scale intervals
cover the full bracket and retain the zero-distance degeneracy.

Retain every interval's bounds, optimization gaps, anchors, split history,
final ratio and stopping status. Compare every result with the old known-noise
and separately profiled-noise bounds. In exact arithmetic, an optimized
common-scale null contains sigma=1, so its optimal log-e-value cannot exceed
the known-noise optimized value; account for the numerical optimization gaps
when checking this. A 30-minute total wall budget is checked between cases.
Any incomplete run remains explicitly incomplete. No scientific power or
real-data calibration conclusion follows just from passing the algebra tests.
