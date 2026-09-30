# Response to focused cubic weight-design audit

This is a technical response, not a full ICML review or acceptance assessment.
The actual `claude-fable-5-1` identity and successful terminal status are in the
unaltered response stream and manifest. The empirical optimization is still
running. R1/R2/R7 from the first full review remain unresolved.

- **W1, numerical operator error:** agreed. The 256-epsilon scalar guard is
  roundoff-only and does not enclose FINUFFT approximation. A declared final-weight
  check compares the nominal Gram at tolerances 1e-12 and 1e-14, records quadratic
  and action differences, and shows heuristic sensitivity margins. It also
  compares one fixed cubic forward/adjoint pair. The runner waits for completed
  weights; no result is available yet. Two approximations can share an error,
  so observed differences will not be called a validated operator bound or used
  to claim an extra probability guarantee. The original fit/audit remains intact.
- **W2, full-objective packing/supports:** added a two-particle, two-frequency
  problem with unequal geometry, angles, shifts and transfers. Its eight real
  coordinates are assembled explicitly; the sinc Gram and target integrals are
  independently computed. The conic operator and remainder components still use
  the tested implementation; this is not an independent physics derivation.
  The new SCS optimum is 1.9982265820361, explicit-gradient value 1.9982265820365,
  dual lower 1.9982256752815 and upper 1.9982265820368 (gap 4.5378e-7).
  Random, non-Ritz unit spatial directions also satisfy full multi-particle
  gradient and unsmoothed norm-support tests. The first new test had a broadcasting
  error in its independent integral code; the second attempt failed in CLARABEL.
  Both logs and exact test versions are preserved. SCS solves the unchanged
  geometry, objective and tolerances declared for its own run; no empirical
  target or case is altered.
- **W3, zero boundary:** a fresh rerun retains `optimal_inaccurate` and records
  2.5289681687701 versus analytic optimum 2.5289681687657, with weight norm
  4.77e-12. This characterizes a rerun; it cannot recover the discarded numeric
  value of the original attempt. The warning remains visible in the test log.
- **W4, selected guide:** the pending final-weight check will report the selected
  squared Ritz value, fresh Rayleigh lower and upper. Best-guide selection may
  miss a candidate with a better certified upper. No candidate is reselected
  and no extra probes are drawn in response to this diagnostic.
- **W5, gap:** the lower is a valid support for the continuous residual objective
  and hence its larger nominal-residual majorant. The audited gap can include
  majorization, entropy, numerical safeguards and spectral-upper looseness; it
  is not a measurement of optimizer error alone. The diagnostic will show eta,
  entropy and roundoff terms without pretending they decompose the whole gap.
  No sum-width gap is transferred to a folded-normal or joint-refined interval.
- **W6, default seed:** the historical implementation is now included in the
  twelve-estimator source lock and remains byte-identical for replay. A new
  `uq_cubic_design_api.fit_with_fresh_certificate` entry point requires the seed
  explicitly and rejects a shared design/certificate seed. Future experiments
  should use this entry point plus their recorded seed inventory. This is
  hygiene, not a replacement for design/probe independence. The current run's
  prospective, uniquely inventoried seed is unchanged; no seeds are hidden in
  NPZ metadata in this run.
- **W7, failure paths:** new checks force an empty second eigensolver return and
  verify that the first emitted history/checkpoint survives. A nonfinite spatial
  mode is rejected before evaluation/certification. The running study continues
  to save errors and all evaluation logs; it does not infer validity from solver
  success or convergence from a valid final audit.
- **W8:** no change to the separate cross-term post-audit or its order-64 pad.

All five added tests pass; their numerical checks do not establish usefulness,
experimental nuisance calibration, or floating-point certification. Final-weight
numerics and the empirical result must still be incorporated after completion.
