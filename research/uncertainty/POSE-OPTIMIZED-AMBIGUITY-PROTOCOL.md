# Prospective pose-optimized constant-cell ambiguity probe

This development protocol is specified before its first outcomes. It extends the
completed selected-pose two-point lower bound, not a confidence estimator or an
experimental coverage measurement. It does not replace the original grid.

Start with a completed v3 fixed-pair witness on the same 128-particle/radius-5
geometry, known white noise, B=2, unit-norm pilot, 0.07-field target, and joint
rotation/translation ball. Project its two unscaled continuous densities onto
the exact orthonormal 24-cell basis. The pilot already belongs to that basis;
orthogonal projection cannot increase distance from it. Keep the original
residual amplitudes, then record the projected pair's class distances and target
gap. Numerical feasibility is guarded inward. This is a subclass of the full
continuous density class, so a feasible testing pair remains a lower bound for
that class. It is not a new assumption on every unknown density.

Hold the two projected densities fixed and minimize their whitened image-mean
squared distance over both collections of allowed poses. Each five-vector is
projected into the unit joint ball. Use analytic constant-cell Fourier first
moments for derivatives and autodifferentiation only through the small rotation
and translation maps. L-BFGS uses at most 80 function evaluations (line-search
checks can exceed this). Its best feasible iterate is retained; local search
need not reach a global minimum. No observed noisy image values are used.

Independently evaluate the saved feasible pair using the existing cell Fourier
forward operator. As in the original construction, scale BOTH full densities
toward zero by min(1,tau/distance_upper), with tau=0.9999*2*Phi^-1(.95).
Because 0 is in the original pilot-centered ball, this keeps both densities in
the class. The scaled feature gap/2 lower-bounds the deterministic half-width
of any uniformly 95% interval under the stated known white-Gaussian law.
The class/pose/target/noise normalization must match the original upper audit.
The continuous target norm remains the no-data denominator, not its cell
projection. Numerical checks do not establish interval-validated arithmetic.

The first pilot is 10049/center/2 degrees/random_boundary. If the implementation
passes independent gradient/forward checks, run the identical procedure on its
saved pair, retain success or failure, and assess whether the larger lower bound
materially changes the interpretation before launching a full grid. No selection
of a favorable result can establish global sharpness or experimental calibration.
A future alternating density/pose optimization would require a separate protocol;
this probe only optimizes poses for a fixed pair of projected densities.

## Pilot outcome and prospective grid-resolution check

The 24-cell pilot completes in 8.59 seconds. Local pose fitting improves its
projected-pair lower width from 0.08362 to 0.09624 of no data, but this is below
the original continuous fixed-pair witness's 0.12634. The projection loses enough
feature energy and changes the predicted means enough to erase the benefit.
This is not evidence that optimizing poses strengthens the existing lower bound.
The optimizer stops at the function-evaluation limit, not a global certificate.

Before a larger pose grid, test nested 48- and 72-cell representations on the
same pilot case and the same 80-evaluation budget. The original 24-cell pilot
function is prolonged exactly, with norm-preserving coefficient scaling; do not
refit or renormalize a new pilot and silently change the class. A regression
checks exact Fourier projections and norm preservation. Both resolution outcomes
will be retained, including negative ones. Record the final densities' negative-
value energy fractions to make the signed-class nature of the witness explicit.

## Resolution outcomes and prospective continuous-density refit

The 48/72-cell probes finish in 80.30/192.53 seconds, with lower widths
0.12169/0.12438; both remain below the original continuous 0.12634 witness.
Their optimized image distances fall below the testing budget, so common scaling
is one. This exposes the fixed-density limitation: once scaling reaches one,
reducing image distance cannot increase that pair's feature gap. Negative-value
energy fractions are about 0.49/0.36 for the two densities, consistent with the
explicit signed class; no biological admissibility is implied.

A new, bounded alternating probe is therefore specified before its results.
Hold the optimized 72-cell poses fixed, return to the FULL continuous-density
modulus solver, and refit its two densities using the original target, pilot,
noise, B/P and testing threshold. Warm-start from the original continuous
weights, use the same order-40 Gram/error pads, rank-1024 preconditioner and at
most 25 outer iterations. Record the constructive lower and fixed-pair dual
upper regardless of convergence. This can use the newly available image-distance
slack; the full continuous solver does not inherit a cell restriction.

Allow at most three density/pose alternations for this single pilot case before
assessing value. Each pose step uses the nested 72-cell projection and 80 function
evaluations; each following density step returns to the continuous class. Preserve
all intermediate feasible witnesses. Keep the best lower, but do not call it a
global modulus optimum. Selection among deterministic feasible testing witnesses
introduces no sampling/coverage event; numerical validity still needs the stated
checks. No new spectral random certificate is drawn or selected by this procedure.
