# Pose-penalty coordinate metric: development only

30 September 2026. The ongoing cubic run preserves its declared nominal-Gram
preconditioner. Its early line search explores large rejected guide values
before a small improvement. This motivates an independently tested coordinate
metric, not changing that empirical run or announcing a new empirical outcome.

For a local sum of norms, drop negative rank-one Hessian terms to form a
positive block metric from the noise norm and quartic Sobolev norms. Add the
nominal density Gram's existing low-rank factor scaled by B/H. The shared
spectral and pilot curvature are omitted: this is not a Hessian upper bound for
the whole objective. Fixed positive floors are merely design choices.

Let the resulting metric be H=D+UU^T, with strictly positive frequency blocks D.
Set K=D^(-1/2)U, R=I+KK^T and T=D^(-1/2)R^(-1/2). Then
T^T H T=I. Low-rank eigendecomposition applies R powers without a dense global
matrix. The forward map is w=T x; the gradient map is T^T grad_w, in the reverse
order. T is generally not symmetric. Initial coordinates use
T^(-1)w=R^(1/2)D^(1/2)w. These are elementary linear-algebra identities, not a
statistical contribution or new coverage guarantee.

The implementation is outside the frozen source lock. No empirical fit uses it.
A new bounded protocol would be required before any further empirical run, and
all current outcomes must remain preserved. Tests compare the whitening,
inverse, chain rule and local metric against independently assembled dense
matrices, including multiple particles, real/imag packing and zero-rank factors.

The implemented spectral factorization discards numerical-null low-rank
directions and reports each cutoff and the largest omitted positive eigenvalue.
The exact whitening identity applies to the untruncated real-arithmetic formula;
finite-precision truncation yields approximate whitening of the supplied metric.
Remainder norms use the larger of the relative median floor and 1e-30 to handle
zero derivative Grams without overflow. Neither approximation affects interval
validity because this component only changes optimization coordinates.

## Bounded compute-only probe, before any further fitting

At the original fixed-pose EMPIAR-10049 pilot_region_1 weights (not any current
optimized checkpoint), build this metric once with rank 1024, one degree/0.5 A,
B=2/P=1 and the saved fixed-pose residual. Record setup time and peak memory.
Use four normalized Gaussian diagnostic pairs with seed 650171 to check inverse
round trips, T^T H T against an independently applied block/low-rank metric,
and the coordinate-gradient duality identity. Record errors and the declared
1e-5 relative threshold; retain failures. This is an operator computation test,
not a fit, a new interval, a coverage check or evidence of faster convergence.
Commit its source before running it. No optimized weights are produced or
substituted into the ongoing protocol.
