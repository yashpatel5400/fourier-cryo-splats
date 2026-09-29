# Frozen continuous-density/pose validation v1

Date: 29 September 2026. Freeze before generating any study signal, pose or noise
draw, fitting any study weights, or inspecting any study coverage outcome.
A separate smoke run used old development geometry, eight particles, different
seeds and only 100 noise draws. It validates execution, not the study hypothesis.

The additional-exposure prediction experiment has already been evaluated. Its
published-pose geometries are the acquisition population here; these are not
new molecular structures or uninspected raw experimental cohorts. This study
provides fresh simulated designs/signals/noise under a locked method. Published
consensus poses and CTFs remain conditional. Unknown experimental density is
never used as a coverage label.

## Fixed design and signal class

Use EMPIAR-10028/10049/10076, four geometry replicates each, 128 particles per
replicate, Fourier half-plane radius five (40 complex pairs). Permute each
additional 4096-image pool with seed 609451 + integer accession; use consecutive
disjoint 128-particle blocks for replicates zero through three. Exposure groups
may recur across replicates; do not portray them as twelve independently drawn
biological datasets. Record exact geometry indices/hashes.

Use the existing unit-L2 24-cell Gaussian-derived pilot from development, with
its checkpoint locked by hash. On the full unit cube H=L2([-1/2,1/2]^3), prescribe
||rho-rho0|| <= 2. Fix the Gaussian whitening scale by the existing pilot signal
power convention (nominal signal-to-noise power 0.1). It is known simulation
noise, not an estimated experimental calibration.

For each dataset/replicate use four exact 64-cell generators, each of unit L2
norm: the deposited-map reference and three newly seeded procedural densities
(unions of ellipsoids, a modulated shell, and a signed multiscale field). The
triangle inequality places all within the prescribed radius two. The procedural
generators are structural controls, not additional measured biological samples.
The scripts fix their distributions and seed formula 609461 + accession*100 +
replicate*10 + generator index. Retain their hashes and the signed control even
if it is unfavorable or less physically representative.

Evaluate central Gaussian averages and axial differences centered at z=+/-0.08,
with standard deviations 0.03 and 0.07 of the physical field. Targets, coordinate
frame and scales are fixed before outcomes. These are pointwise feature tests;
no simultaneous confidence-map or data-selected-target claim is made.

## Locked procedure and controls

For each target, optimize the continuous fixed-pose sum-width objective at B=2
using Gauss order 40 and rank-1024 pivoted-Cholesky preconditioning, maximum 100
outer steps and 0.5% objective-gap tolerance. Preserve nonconverged cases and
report their still-valid upper bounds; never discard by convergence or width.
Independently recheck the continuous residual with analytic sinc integrals.

Post-audit the same weights for rotation budgets 0.1 and 0.5 degrees with joint
five-dimensional pose unit balls and shift radius 0.01/24 field units. Use the
continuous polynomial-Fourier spectral audit, Gauss order 32, explicit quadrature
padding and cubic remainder. Do not claim these weights optimize nonlinear
width. Compare five procedures on each identical signal: noise-only interval,
24-cell fixed-pose audit, continuous fixed-pose audit, continuous nonlinear-pose
audit, and the no-data bound. Each chooses the no-data center/width if its nominal
width is larger, using design only. Mathematical covariance and group-bootstrap
comparators remain separate completed development studies, not unperformed
external code reproductions.

For each of the four generators, evaluate nominal poses, independent random
pose-ball boundary points, a coherent x-rotation boundary, coherent pose radius
exceeded by factor four, and Gaussian noise scale exceeded by factor two. The
last two deliberately violate assumptions and must remain separate from covered
cases. Pose draws use 609481 + accession*100 + replicate*10 + 10*angle_degrees.

Also include both signs of the exact full-L2 nominal bias boundary and both
signs of exact density-ball elimination at the coherent nonzero pose. The latter
uses independent sinc integrals and exact constant-cell pilot projections; it
is a feasible lower bound on the global pose supremum, not a global adversary.
These four boundary cases are inside the declared model. No generator is
restricted to the inference quadrature nodes by assertion.

## Outcomes and analysis

There are 96 audit settings and 2304 signal/scenario records: 1536 in-model and
768 deliberate out-of-model controls. Report exact Gaussian noise-marginal
coverage, relative half-widths, no-data fallback counts, sign-certification
probabilities, optimization gaps, feasible-bias/upper-bound ratios, timing and
stored-array memory. Summarize all targets/radii/geometries without choosing a
favorable subset. Separate typical generator controls from estimator-specific
bias boundaries. The main correctness criterion is minimum in-model analytic
coverage >=0.95 within numerical tolerance 1e-8. Correctness alone does not
establish practical usefulness, novelty or experimentally calibrated assumptions.

Cross-check each record using 10000 scalar independent Gaussian noise draws
with stream seed 609471 + accession*100 + replicate. Pair draws across methods
within a case and report marginal exact-binomial 95% intervals. These draws are
replications of measurement noise, not independent molecules, voxels or datasets.
No multiplicity-adjusted map-level coverage claim follows from these checks.

## Integrity

The lock manifest records exact code, pilot, reference and acquisition-metadata
hashes. Commit it before execution. Archive executed source bytes and environment.
Fail on hash drift and preserve all completed cases. Numerical fixes require a
dated amendment and retained error/output; scientific changes after outcomes
belong to a separately labeled study. There is no outcome-dependent stopping.
All three datasets and four geometry replicates must finish, including failures
of usefulness. This study is one component of an ongoing research candidate;
independent scientific review and broader practical validation remain required.
