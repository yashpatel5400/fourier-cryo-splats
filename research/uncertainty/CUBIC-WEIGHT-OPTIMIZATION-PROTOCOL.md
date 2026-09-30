# Bounded pose-aware cubic weight-design probe

Development protocol, 30 September 2026, after the completed cubic fixed-weight
case. That case's bias is 2.46666: 0.81489 remainder and 1.65177 polynomial/density
terms, with noise adding 0.60113 to half-width. Sign power is 1.48e-10. This
motivates changing weights, not another remainder-only claim. This is not frozen
confirmation; the case was selected during method development.

## Single empirical case and fixed inputs

Use only EMPIAR-10049 pilot_region_1, 128 particles, radius 12, the locked 20 A
Gaussian target, B=2/P=1, one degree and 0.5 A in each particle's joint pose ball.
Start from the completed fixed-pose weights. Keep the original continuous pilot,
noise convention, indices and all original reference scenarios. Reuse the
positive degree/particle scales from the completed cubic audit; these are design
inputs, not its random certificate. No new particles, reference target choice,
pose radius, density class or target-width selection is part of this probe.

## Objective and bounded algorithm

With fixed scales and lifted radius L, optimize the convex triangle surrogate

  z ||w|| + B ||ell-A0*w||
    + B L sqrt(||F_Q(w)||_op^2 + E6 ||M a(w)||^2)
    + L ||F_cont(w)*rho0|| + (B+P) r_PSD(w).

Here F_Q is the full cubic field at quadrature order 64, E6 is its degree-six
analytic kernel remainder, M maps real/imaginary pair amplitudes to absolute
scaled column coefficient sums, and a(w) is that amplitude vector. The known
pilot pairing uses twenty analytic continuous cell moments. r_PSD is the convex
Bell-weighted sum of padded real derivative-Gram norms plus the approximate
embedding residual pad, as implemented in uq_sobolev_penalty.py. The fixed-pose
residual uses order-80 continuous quadrature. z is the two-sided Gaussian
critical value at alpha_noise=.05/12-1e-6/12. This sum-width surrogate is not
asserted to be the exact folded-normal optimum.

Guide L-BFGS-B with three approximate leading Ritz modes (tolerance 1e-3,
subspace size 13, at most 100 eigensolver iterations), log-sum-exp smoothing
0.002, and warm starts. Partial returned modes may guide design but never
certify an upper bound; no returned modes is a retained numerical failure.
Use at most 30 optimizer iterations, requested maximum 40 function evaluations,
maxls=10, ftol=1e-5 and gtol=1e-6. The library checks its evaluation budget after
a line search, so record any overshoot up to that line-search allowance. Use
optimization seed 650011. Preserve the lowest observed surrogate candidate,
all evaluations, partial-mode states, checkpoints, and optimizer status.
Use the existing continuous Gram preconditioner with rank 1024, with its initial
ridge fixed by the saved starting weights; report setup time separately.

Final certification uses the SAME order-64 pose surrogate, with four fresh
Gaussian probes, forty steps, seed 650001 and delta=1e-6/12. Verify before the
run that this seed has no prior certificate use. The objective's lower bound,
if reported, must use valid continuous density supports and pose/remainder norm
supports for this same order-64 surrogate. A Ritz value or optimizer-success
flag is not a certificate. Do not compare a lower bound for order 64 with an
upper bound for a different surrogate at order 80.

After certification, the existing joint residual/pose cross-term bound and
known-pilot refinement may tighten the interval at the same weights/event.
An independent direct Bell/Sobolev evaluation may provide a smaller deterministic
remainder. Keep the surrogate gap distinct from this post-audit width. No-data
fallback changes both center and width. The old quadratic, cubic fixed-weight,
and new optimized intervals are separate alternatives; no unadjusted selection
minimum is claimed as a joint confidence interval.

## Prerequisites and outputs

Before empirical fitting, commit source and this protocol; test full objective
weight gradients and supporting inequalities on independent toy geometry, and
compare a small problem with a separately assembled conic spectral-norm solve.
Preserve any failed test or numerical attempt. The new operator/penalty component
tests are implementation checks, not experimental calibration.

Persist exact source/input hashes, geometry/weight arrays, positive scales,
integration pads, spectral history, error allocation, original/refined bias
components, width, sign power and coverage for every original reference scenario,
runtime, memory, optimizer state and a valid surrogate gap if available. Stop
this bounded study after its one final outcome, favorable or not. Additional
cases, iterations or starts require another explicit development declaration;
they must not overwrite or be represented as this probe's result.
