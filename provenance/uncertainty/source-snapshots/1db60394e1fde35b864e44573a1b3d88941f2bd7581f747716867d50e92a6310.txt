# Local Gaussian pose-marginalization comparison

Development protocol, 30 September 2026, before these outcomes. The completed
Fourier baseline supplies useful broader-prior comparisons but holds poses
fixed. This extension checks one missing approximation: local Gaussian pose
uncertainty. It does not solve the full review's experimental-calibration problem.

Rangan et al. (2024), methods summary and equations 2–12, analyze coupled
volume/alignment curvature and account for optimal alignments changing with
the volume. Their soft-mode calculation removes rigid rotational directions.
Our comparison below is a simpler pilot-linearized Gaussian model, not their
volumetric Hessian algorithm, its complete nonlinear curvature, or a reproduction
of their experiments. Primary: https://arxiv.org/html/2411.13263v2.
The existing fixed-pose Fourier baseline independently specializes Ullrich
et al.'s Gaussian objective, with the documented project extensions:
https://proceedings.mlr.press/v115/ullrich20a.html.

## Exact statement for the approximate model

With the existing finite Hermitian acquisition A and pilot coefficients m,
use y=A x+J u+epsilon, x~N(m,tau^2 I), u~N(0,I), epsilon~N(0,I), independent.
Each particle's five columns of J are the derivative of the *same trilinear
pilot signal*, scaled by the declared coordinate standard deviations. Rotation
uses k exp([omega]_cross); detector translation multiplies by exp(-2 pi i q.t).
Reject non-differentiable interpolation knots instead of hiding a one-sided
derivative. This linearization omits density/pose products and higher orders.

Marginal noise covariance K=I+J J^T has independent particle blocks. Whitening
with W=K^(-1/2) gives Q=A^T K^(-1) A+tau^(-2) I. For a feature ell,
v=Q^(-1)ell and w=K^(-1)A v give posterior variance ell^T v. The identity

  ell^T v = ||w||^2 + ||J^T w||^2 + tau^2 ||A^T w-ell||^2

is ordinary Gaussian conditioning, not a new theorem or a uniform coverage
result. It follows from Qv=ell and A^T w-ell=-tau^(-2)v. The implementation
uses small Gram eigenproblems for W and the existing matrix-free Gaussian
solver. Numerical tests compare against a separately assembled dense joint
density/pose posterior, and finite differences of the nonlinear interpolator.

## Bounded empirical comparison

Use all twelve locked features, the same 128 particles, radius 12, geometry
seed 609315, 33^3 Fourier grid, pilot, known-noise convention and independent
continuous reference generators as the earlier matched baseline. Retain both
coordinate prior SDs 0.1 and 1.0; these are the previously inspected favorable
broader priors, so this is post-outcome development, not a new prior selection.
For each prior, compare fixed poses with per-coordinate Gaussian SDs of one
degree for each rotation coordinate and 0.5 A for each detector shift.
Gaussian norms are unbounded; rotation RMS is sqrt(3) degrees, not a one-degree
hard joint pose ball. The external target frame remains fixed by the prior;
do not project away a gauge mode while keeping an unchanged spatial target.

This gives 48 target/prior/pose-model solves across three geometries. Each uses
CG rtol 1e-10 and maxiter 2000. Record convergence failures and stop that dataset
without overwriting; attempt the other datasets independently. Recompute the
fixed-pose controls and compare with the exact saved weak-prior outcomes.
Keep alpha=.05/12 per feature for each fixed prior/pose model; no guarantee is
claimed over selecting hyperparameters. Record posterior widths, raw measurement
SD, prior variance identity, and the eight original reference scenarios with
their analytic fixed-signal coverage and sign power. This is conditional known-
noise verification and sensitivity, not coverage measured on unknown density.

At the same prescribed coherent/random 0/1/2-degree boundary configurations,
also report the pilot's nonlinear-trilinear versus linearized forward discrepancy.
This exposes local approximation error separately from the previously measured
continuous-reference interpolation error. No perturbation or target is selected
by these outcomes. Save raw weights, Jacobian, exact hashes, timing and memory.

Commit source and this protocol before fitting. Preserve every outcome under
`pilot-selected-fourier-pose-v1`. No further priors, pose scales or restarts are
part of this comparison. The supplied pose prior, noise covariance, homogeneity
and finite representation remain assumptions; this is not an ab-initio run or
an independently calibrated pose posterior for the experimental particles.
