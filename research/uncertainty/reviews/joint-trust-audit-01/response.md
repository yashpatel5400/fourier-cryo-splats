# Response to focused joint trust-region audit

30 September 2026. This is a response to the authentic mathematical audit,
not a new ICML acceptance assessment. The first full review remains reject.
The queued empirical job was stopped during its scheduling wait, before any
result directory or joint empirical fit existed. The cancellation is retained.
The separately running adaptive triangle-objective fit was not changed.

- **T1:** Each failed shifted Galerkin solve is now recorded. A failed scalar
  evaluation returns the valid old-cross value to the scalar search, and the
  final bound keeps the minimum over successful resolvent values and that old
  value. If all solves fail, `selected` is null and the old bound remains.
  A singular-cache test initially assumed LAPACK would fail at every shift;
  some floating-point factorizations returned values. That test failure is kept.
  The corrected test checks this mixed case and separately injects solve
  failures to exercise the actual all-failed branch. Returning a finite valid
  scalar fallback also avoids the earlier scalar interpolator's infinity warnings.
- **T2:** Every evaluated coefficient vector is rescored against the full final
  cut set, retaining its own oracle guide. The selected gap is measured against
  the last master, and the last iteration gap is separately labeled. The selected
  final weights are checkpointed even if they were not an earlier improvement.
  A test makes an optimistic early candidate lose after final-cut rescoring.
  The resulting scores remain approximate, padded quadrature guides; an unseen
  adverse pose direction can still be missed. No full-space convergence is claimed.
- **T3:** The secular equation is solved in log shift with stable log-sum-exp
  norm evaluation. Tests check the norm before boundary extension at three
  near-hard scales; the feasible boundary extension remains. This fixes root
  tolerance without changing the ball or any empirical parameter.
- **T4:** A three-column/two-start Krylov test with U=1.5 lambda_max gives
  upper 16.17582 versus old-cross 23.11710, retaining the lower exact witness
  check. This small case does not predict production improvement. A pure joint
  test zeros pilot and remainder slack and compares to order-16 direct
  integration: reference 12.00744215539, bound 12.00744215678. The independent
  robust-norm SDP now explicitly establishes independence only for the joint
  term; pilot/remainder matrices are shared. Accepted Krylov origins/depths are
  recorded, so two-start column counts cannot be misread as per-start depths.
- **T5:** Keep the conservative PSD factor guard and report padded and unpadded
  cut norms separately. Removing the target diagonal pad could undermine the
  generic PSD correction; it is not necessary for this correction.
- **T6:** Both public routines validate alpha/delta, the corresponding critical
  value and distinct nonnegative integer design/audit seeds. A separate upper-
  component helper avoids the unused triangle-dual calculation, while frozen
  older source stays unchanged. A regression replays every shared upper field
  against the prior implementation.
- **T7:** The theory note now explicitly retains the conservative degree-six
  cross pad and reduced amplitude-triangle remainder. Neither was loosened
  to obtain favorable results.

The corrected full suite passes 246 tests, with one skip and the earlier
intentional conic warning, in 184.09 seconds before the first empirical joint fit. The original successful/failed prerequisites and original review
packet are preserved unchanged. This response does not claim that Fable has
reviewed these corrections or endorsed publication readiness.
