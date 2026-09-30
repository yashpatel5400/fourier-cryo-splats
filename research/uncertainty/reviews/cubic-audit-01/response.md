# Response to the focused cubic mathematics audit

This is a response to a verified `claude-fable-5-1` invocation, not a full ICML
review. The first full review remains a rejection. The empirical cubic audit
was still running when these checks were performed; its weights, scales, seed,
iterations and error allocation have not changed.

- **C1:** checked all 28,160 actual frequencies and all twenty Gaussian moments
  against independent separable 200-node quadrature at the actual field size
  236.1600037 A, sigma 0.08468835 field units and locked center. The largest
  per-moment-column relative norm error is 2.346e-14, below 1e-10. This is a
  numerical check in the audit regime, not a global special-function guarantee.
- **C2:** the exact numeric seed inventory covers 515 result JSON files and finds
  640001 only in the current cubic case. The source quadratic seed is 620101.
  A whole-word search of pre-declaration commit 9a9725f finds no 640001 in
  source, protocols or results. The initial verification attempt falsely failed
  because a substring search matched digits inside unrelated decimal results;
  both that failed record and its exact archived code remain available. The v2
  check changes the search to whole words. No certificate was redrawn or chosen
  after seeing an outcome.
- **C3:** the estimator is explicitly clarified below and in the manuscript.
  Independently reconstructing the continuous fixed-pose residual, its Gauss
  integration pad and the original arithmetic guard gives bias
  0.0121459269552 versus stored 0.0121459269564 (relative error 9.633e-11).
  The exact original active runner is preserved rather than silently changing
  it mid-run; the completed regression is a separate prerequisite for interpreting
  its output.
- **C4:** the requested mixed planar test was added. A proposed factor-20 slack
  threshold failed: the ratio is 35.57084 with exact planar geometry and with
  1e-7 distortion. This is retained as a conservatism finding, not concealed by
  labeling a relaxed threshold as success. The mixed test now checks validity
  and independent fourth-order error scaling; a separate pure-translation test
  has a bound/error ratio below three and the expected factor-sixteen scaling.
  The scalar nonlinear rotation recurrence separately detects an omitted Phi3.
  These checks strengthen falsification, but do not establish a tight general
  rotational remainder. Both the failed and subsequent test logs are retained.
- **C5:** an additional analytic pair-frequency calculation evaluates the
  continuous polynomial-field block traces for the fixed particle positions
  0, 64 and 127. Across their nine degree blocks, continuous norm / design
  scale is 0.9823--1.0285. This selected diagnostic does not indicate a large
  aliasing distortion, but does not prove optimal scales or inspect every
  particle. The exact record is `audit-regressions/cubic-design-scale-check.json`.
  It uses continuous polynomial cube moments independently of order-12 design
  quadrature. No scale retuning is applied to the active random event.
- **C6:** one independent actual order-80 operator adjoint check has relative
  pairing residual 5.98e-16. This cannot upper-bound the operator approximation
  error, so it is reported as a heuristic check rather than misused as a rigorous
  spectral margin. The lack of validated arithmetic remains explicit.
- **C7:** the follow-up record and manuscript explicitly label the quadratic
  and cubic intervals as alternatives. Their unadjusted minimum is not claimed
  to be one joint confidence interval. The locked twelve-feature family still
  uses its originally declared quadratic procedure.

Completed numerical evidence is
`results/uncertainty/development/audit-regressions/cubic-audit-01-checks-v2.json`.
The earlier failed checker is `cubic-audit-01-checks.json`. The cubic audit subsequently completed at 2,187 seconds and 1.42 GB peak
resident memory. Relative width is 0.2523 versus the quadratic enclosure's
0.4734, but minimum reference sign power is only 1.48e-10. Its bias bound is
2.4667, including 0.8149 remainder and 1.6518 polynomial/density terms.
These are later scientific results, not facts the focused reviewer had read. The polynomial trace helper was
checked against independent tensor quadrature. One initial scalar test tolerance
was tighter than quadrature cancellation at a known exact zero (5.8e-15); the
failed log is retained and the corrected test separately checks that exact zero.

## Clarification of the original conditional proposition

Use the pilot-corrected estimator

    theta_hat = <ell,rho0> + w^T (y - A0 rho0).

For rho=rho0+g, the noiseless error is

    <A0* w-ell,g> + <(Au*-A0*)w,rho0+g>.

Thus no uncorrected term <ell-A0*w,rho0> belongs in its bias. This clarifies
the original estimator used by the code and does not change the protocol,
statistical class or reference centers. The stored fixed-pose bias is B times
the continuous residual norm with the solver's recorded integration/arithmetic
padding; the completed independent reconstruction checks that interpretation.
