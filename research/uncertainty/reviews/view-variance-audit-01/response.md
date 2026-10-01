# Response to focused view-variance audit 1

The authentic requested `claude-fable-5-1` audit is preserved unchanged. It was
prepared from source 8c324d8, checked the supplied text/code by reading, and
executed no tests or numerical arrays. Its packet was unchanged during the
invocation. This is not a fourth full review. All three full reviews still
recommend rejection. The new analyses below were motivated by this audit.

- **V1, fixed assumption statement:** the new [comparator protocol](../../paired-power-v1/VIEW-RISK-BASELINES-PROTOCOL.md)
  and paper explicitly require physical amplitude conditionally independent of
  noise given view for grouped methods. They distinguish the larger null
  protected by the individual envelope. A finite, exactly enumerated cubic
  counterexample has grouped mean .75 versus noise-adaptive event probability
  one; every fixed amplitude has probability at most .5. This corrects the
  omitted qualifier without claiming its experimental validity.
- **V2, retained:** eta_L is explicitly the finite-replica conditional envelope
  mean. The shared confidence event and binomial argument are unchanged. The
  chi-squared-only formulation discards the density-cap branch; no experimental
  radius or ease of estimating it is asserted.
- **V3, implemented:** two classical CVaR confidence comparators (uniform DKW
  and independent split empirical Bernstein), a mean-only variance envelope
  and an unpaired variance envelope are complete on all eight candidate/score
  cases. All 480 projections are retained in [results](../../paired-power-v1/VIEW-RISK-BASELINES-RESULTS.md).
  A separate linear program checks the CVaR optimum, including fractional
  weights at boundary atoms. The paired method is tighter at kappa=1.1 for
  this saved budget, but this is not a population CVaR improvement or proof
  against better simulation allocation.
- **V4, confirmed:** all raw product estimates and concentration components
  are now tabulated. The removed 10049 raw estimates are negative, while
  concentration slack produces positive variance upper bounds. We no longer
  interpret these bounds as measured population viewing sensitivity.
- **V5, partially addressed:** the paper distinguishes projection/noise budget
  matching from the older Clopper--Pearson allocation. Paired critical-count
  differences are retained for every comparator. Calibration repetition and
  an allocation frontier remain unperformed. Single-image pointwise intervals
  are not a paired test of method ordering or calibration variability.
- **V6, additional controls:** the [fresh preferred-view study](../../paired-power-v1/PREFERRED-VIEW-RESULTS.md)
  evaluates all three axes, kappa=1.1/2/5 and amplitudes .9/1/1.1, retaining
  both correct-null maps and every alternative. It records actual 64 disjoint
  groups of 1,024 images as well as all projections. Sensitivity changes
  substantially with amplitude. This is not unconditional type-I validation
  over calibration repetitions. The original report's printed zero can mean
  exact non-rejection or numerical tail underflow; the full critical counts
  distinguish these, and new nonzero-tail summaries avoid printing underflow
  as an exact probability. The audit did not see n=100,000: at kappa=1 the
  10028 combined projection there is .782 [.175,.994], so its broad opening
  statement of no power on 10028 must be restricted to the displayed n=10,000
  setting. This does not establish reliable power.
- **V7, tests added:** exact-edge, nearly flat stationary and tiny-cubic cases
  pass, along with the assumption counterexample. These are ordinary
  floating-point checks, not certified interval arithmetic. Near-zero
  discriminants remain a rounding issue; a mathematical double root should
  not be used to declare all rounded near-double roots harmless. A repeated
  calibration coverage experiment remains outstanding.
- **V8, broadened verification:** a separately seeded 1% sample (328 views per
  stack) reproduces 167,936 cubics and 10,496 grouped/individual counts using
  separable physical-cell sums, independently coded features, Vandermonde
  coefficients and companion roots. This is broader than the original first-
  view replay, not a proof covering all views or rounding. The preferred-view
  verifier explicitly checks unique nonzero frequencies with no antipodal
  pairs, and exact equality with the first archived transfer profile. It
  replays every new event count, critical value and repeated-group vector.
- **V9, unresolved practical calibration:** the new cap laws are synthetic and
  supplied exactly. No known bound relative to Haar or an estimated viewing
  law is established for experimental particles. An estimated-Q construction
  would still require credible divergence calibration. Larger ratios remain
  weak in these controls. The suggested event-probability extremal tilt has
  not been run; cap laws are not represented as that experiment.

The most consequential open issues remain practical non-oracle scores,
experimental noise/amplitude/view calibration, independent calibration
repetitions, extremal-law controls, smaller structural changes and useful
results across all three stacks. These additions improve the comparison and
falsification record, but do not resolve the full review's novelty/usefulness
objections or warrant declaring the research goal complete.
