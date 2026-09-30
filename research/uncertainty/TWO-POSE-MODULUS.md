# Prospective two-pose ambiguity test

This protocol precedes its new solver and experimental outcomes. It addresses
whether a wide pose-sensitive interval is computational conservatism or an
intrinsic ambiguity under the declared class. It is a specialization of the
classical two-point testing/modulus argument, not a new general minimax theorem.

## Class and lower bound

Keep the existing continuous cube class `||rho-pilot|| <= B`, with `P=||pilot||`
and `B>=P`, and known white Gaussian observation noise. Select two feasible
pose configurations independently of inference noise, with operators A0,A1.
They need not share a density. For `p=A1 pilot-A0 pilot`, the feature-gap modulus
under observation distance at most tau has convex dual upper objective

    G(w)=B(||ell-A0*w||+||ell-A1*w||)-w'p+tau ||w||.

For any w, let hj=ell-Aj*w and aj=B/upper(||hj||). Construct

    rho1=pilot+a1 h1,  rho0=pilot-a0 h0.

Both are in the class. Their mean difference is
`p+a1(A1 ell-A1A1*w)+a0(A0 ell-A0A0*w)`.
An upward integration-error pad bounds its norm by D. Scale **both densities**
by `lambda=min(1,tau/D)`. This is still in the original class because zero
belongs to the ball (`P<=B`) and the ball is convex. The observation distance
is at most tau, and the feature gap is computable analytically. Every such
pair supplies a feasible modulus lower bound, even if optimization fails.
The dual objective and constructive lower bound bracket the fixed-pair
modulus; optimizing poses is not required for a valid lower bound.

Choose `tau=.9999*2*Phi^{-1}(1-alpha)` with alpha=.05. For whitened Gaussian
means separated by D<=tau, the minimum sum of testing errors is
`2 Phi(-D/2)>2 alpha`. If an interval with deterministic length less than the
feature separation covered both parameter points with probability at least
1-alpha, membership of the first feature would define a test with both errors
at most alpha, a contradiction. Thus half the feature gap is a necessary
half-length for **any** uniformly valid deterministic-length interval on the
larger pose/density class. This does not bound variable-length intervals or
establish experimental nuisance calibration.

## Numerical representation and checks

Use the existing continuous Fourier Gram, analytic Gaussian target and known
constant-cell pilot. Add the existing integration remainder to each residual
norm and to each observation-vector norm. Translation acts by orthogonal
two-coordinate Fourier rotations, so it can wrap a rotated-coordinate Gram
without changing Euclidean error bounds. The unknown densities are continuous
target/adjoint fields plus the pilot; they are not restricted to a voxel grid.
Floating-point/NUFFT arithmetic remains numerically checked, not interval-
validated. Save witness weights, amplitudes, common scale, pose arrays, norm
pads, both objective bounds, all solver statuses and source hashes.

Independent tests will reconstruct witnesses in unrelated small Hilbert
spaces, verify class membership and image distance, compare the modulus
against a direct conic program, and compare the phase wrapper with nonlinear
constant-cell projections. A failed check is recorded before stopping.

## Initial development grid

Use the existing 128-particle/radius-5 geometries, broad .07-field center and
contrast targets, B=2, P=1 and each existing noise scale on all three stacks.
The initial configurations are: identical nominal poses (solver control),
antipodal coherent x rotations, and antipodal random joint-ball boundary poses
at one and two degrees, with translation radius 0.5 Angstrom. Use the fixed
seed 610341 plus accession; the random unit vectors are shared across angular
budgets. Retain every case and compare lower bounds with the pose-aware upper
widths only for matching classes, geometry and nuisance budgets. Since finite
pose candidates do not exhaust the class, a small lower bound cannot prove
the upper bound is loose. This is development after the first full review,
not a new frozen confirmation or a replacement for the failed usefulness tests.

## Numerical hardening after the focused audit

The first complete 30-case run is retained in `two-pose-modulus-v2`. Its source
was reviewed independently by the actual Fable 5.1 model in
`reviews/modulus-audit-01`. The revised implementation rejects negative padded
residual squares and nonfinite witnesses rather than silently clipping them.
It records separate weights and iteration indices for the best constructive
lower and dual upper, the pilot-mean distance, the weight norm, the unnormalized
no-data half-width B||ell||, and a hash of the witness arrays. The Gaussian
testing bound still assumes identity covariance in packed real coordinates.

Ordinary floating-point magnitude pads now move the feature gap downward and
the mean norm and dual objective upward, in addition to the analytic integration
remainders. They are heuristic arithmetic guards, not validated rounding bounds
for the NUFFT, special functions or every accumulation. Exact-continuous tests
with deliberately crude quadrature exercise both norm-pad paths. Additional
tests cover the wrapped target/preconditioner, nonlinear pilot projections and
a case where the best upper and lower occur at different iterations. The same
prespecified grid reruns in `two-pose-modulus-v3`; old outcomes are not replaced.
