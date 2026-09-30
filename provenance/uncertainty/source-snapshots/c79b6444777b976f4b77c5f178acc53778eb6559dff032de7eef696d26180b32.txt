# Conditional bounded coordinate follow-up

30 September 2026, declared while the original cubic weight fit is running.
Early line-search trials reached guide values 928 and 71 before falling below
the initial value 3.045. At the last inspected evaluation 15, the best guide
was 2.626; no final audited interval or gap was available. This motivates a
coordinate change, not a scientific claim about its eventual performance.
The compute-only metric probe and independent two-particle conic check are
complete. This is development, not prospective confirmation on new data.

## Eligibility and one-run limit

Wait for the original `cubic-weight-probe/10049-pilot_region_1-1.json` to finish
its final audit and nonlinear reference checks. Preserve it unchanged. If it
fails numerically, do not silently substitute this procedure. If its reported
relative sum-objective gap is at most 0.005, save a skipped outcome and do not
run another empirical fit. Otherwise execute exactly one new fit under the
settings below. Record the original complete outcome hash and gap. This rule
is declared before the original fit's final outcome is known.

## Identical problem and bounded coordinate change

Use the same EMPIAR-10049 first pilot region, 128 particles, radius 12,
20 A Gaussian target, B=2/P=1, one-degree/0.5 A joint pose balls, original pilot,
noise convention, positive cubic scales and starting fixed-pose weights. Do
not initialize from the original optimized checkpoint. Keep the same order-80
density majorant, order-64 cubic field, analytic pilot moments, quartic
Sobolev penalty and smoothed three-mode design objective.

Replace only the fixed nominal-Gram coordinates by the tested block plus
low-rank metric described in CUBIC-PRECONDITIONER-DEVELOPMENT.md. Evaluate its
density residual at the original starting weights using the same objective.
Use rank 1024, relative norm floor 1e-8 and absolute floor 1e-30. The metric
stays fixed throughout optimization and omits spectral/pilot curvature.
Record setup time, rank truncations and block scales. The transform and its
transpose are different; the new solver uses the latter for gradients.

Preserve the original budget: L-BFGS-B at most 30 iterations, requested maxfun
40 with maxls 10, ftol 1e-5, gtol 1e-6, maxcor 15. Record any library line-search
budget overshoot. Retain every evaluation and every improving checkpoint.
Use three Ritz modes, tolerance 1e-3, ncv 13, maxiter 100 and smoothing 0.002.
Reuse design seed 650011 for the coordinate comparison; this is not a final
certificate stream. Partial Ritz modes can guide design, and empty modes fail.

After choosing the lowest observed guide candidate, compute a new numerical
upper with four fresh Gaussian probes, 40 steps, seed 660001 and delta=1e-6/12.
Check the development JSON seed inventory before fitting. The matrix selection
must never depend on those final probes. Keep alpha=.05/12 and the same
order-64 surrogate on both sides of the primal/dual comparison. Independent
post-refinements and the three original nonlinear reference scenarios are
unchanged. Record all resulting widths, power, gaps, failures, time and memory.

## Interpretation and preservation

This comparison tests whether a fixed coordinate metric improves a bounded
solver on one already selected case. It introduces no new statistical
guarantee, biological target, density class or calibrated pose input. Both
original and follow-up intervals remain separate alternatives; do not select
their smallest width and claim an unadjusted joint confidence level. A small
surrogate gap would resolve optimization of that surrogate, not its bias-bound
tightness or experimental usefulness.

Commit this protocol, runner and all dependencies before execution. Write to
`cubic-coordinate-probe`, retaining the original directory and all failures.
No further starts, longer budgets or changed targets are authorized by this
protocol. The overall research may propose a separately documented next step
based on the resulting evidence, but must not overwrite this comparison.
