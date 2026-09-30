# Response to the focused mathematical audit

This is an authentic Claude Fable 5.1 mathematics audit, not a second full-paper
review or an acceptance decision. The complete stream and unchanged report are
preserved. Round 1's ICML rejection still applies.

- **B1:** The main model and THEORY.md already specified the joint scaled
  five-dimensional ball. That material was accidentally absent from the focused
  packet. The short proof and audit metadata now state the set explicitly,
  including its exclusion of simultaneous maximal rotation and translation.
  A regression test demonstrates that the same lifted radius fails for the
  product set. A separate product-set post-audit expands the first/second
  lifted blocks by their correct factors; it does not relabel older results.
- **B2:** Failed post-audits now save `complete=false`, and summaries reject
  error/numerical-failure records even if a legacy complete flag is true.
  The same status issue in nonlinear diagnostics and the exchange probe was
  fixed. Existing archived cases did not trigger it.
- **B3:** No old-objective optimization gap is paired with a joint width.
  The proof now explicitly distinguishes the objectives and feasible bias
  bracket. The post-audit is not asserted convex or optimized.
- **B4:** Each new post-audit records its source certificate seed and the rule
  that no minimum over repeated certificates is selected. The summary recovers
  older seeds only after verifying the source fit hash.
- **B5:** Summaries list every attempted post-audit, including incomplete or
  failed records. Fallback diagnostics switch to the pilot center with zero
  variance; old files lacking the raw mean are reconstructed when necessary.
  This does not turn unsuccessful optimizations into completed experiments.
- **B6:** New targeted tests compare width-0.07/frequency-5 moments with an
  independent 200-node rule, the order-32 cross vector with order-96 direct
  quadrature against its analytic pad, and a 2-degree cubic remainder with
  its leading rotation direction on the maximizing cube diagonal.
- **B7:** Future fit records save and hash every block scale and denominator,
  and consumers compare the full vectors. Legacy fits cannot retroactively
  gain missing metadata: reconstructed vectors are saved, but this limitation
  is explicitly recorded rather than calling a scalar check a full match.
- **B8:** A new small independent conic fixture uses nonempty spectral cuts and
  an incomplete continuous density projection, checking the full-space lower
  bound against the exact unprojected optimum and global support inequalities.

The targeted regression suite passed 14 tests and the full suite passed 88
tests (`logs/uncertainty/bound-audit-{regression,full}-tests.log`). The separate
product-ball audit of the 10049 two-degree restricted-cut estimator completes:
its sharp-refined relative width rises from 0.5716 in the joint ball to 0.6462
in the product set. This is a larger class, not a loss of validity of the
original joint-ball result. Reference sign power remains negligible. Subsequent
design-scale tests separately check that underresolved quadrature used only to
choose positive scales does not substitute for the final operator evaluation.
