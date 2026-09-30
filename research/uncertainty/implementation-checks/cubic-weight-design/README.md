# Cubic weight-design prerequisites

These are implementation checks before any empirical fit, not evidence of
experimental coverage or power. The first test/code bytes and failed pytest
log are retained. At the original toy noise scale 0.7, CLARABEL returned
`optimal_inaccurate`; the strict status assertion failed. The analytic condition
B||A ell||/||ell|| < z proves that example has the zero-weight global optimum.
A separate regression now explicitly checks that boundary. An informative toy
with the same geometry and noise scale 0.1 has a nonzero optimum and a successful
independent conic solve; the change is a numerical-test choice, not data tuning.

The successful conic value is 2.0715692154662158. The explicit-gradient optimizer
returns 2.0715692154591783; its density-support dual lower is 2.07156919563574 and
audited upper is 2.071569215459424 (relative gap 9.5694e-9). The conic program
assembles the continuous cube Gram directly from sinc kernels and expresses
spectral, pilot-pairing, integration and derivative-Gram terms with independent
CVXPY cones. Common field null spaces are compressed only below 1e-13 relative
singular value, with a reconstruction error assertion. Component field physics
is independently checked in the cubic operator tests.

Directional gradients, convex supporting inequalities (including smoothed
mode supports for the unsmoothed objective), deliberate partial Ritz modes and
the independent final upper audit are also tested. All five checks pass. The
initial, v2 and v3 logs are in `logs/uncertainty/cubic-optimization-prerequisites-*`.
Ordinary floating-point accuracy checks do not certify numerical operator error.
