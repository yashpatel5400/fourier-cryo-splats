# Evidence supplied for full review round 4

This requests a fresh evaluation, not a favorable verdict. Reviews 1-3 remain
unchanged rejections; prior concerns are not declared resolved by new run counts.
The current paper has 33 pages. Its earlier alignment/representation audit
remains largely classical and its practical pose intervals uninformative.

Changes since full review 3:

- The alignment-bias reanalysis reports the estimated-pose failure and compares
  true-pose fits and realized class envelopes; failed refinements remain public.
- Continuous-pose falsification invalidates several moment/covariance candidates;
  all failed, zero-direction and nonconverged cases remain available.
- A separate known-simulator, bounded-view candidate-test branch uses finite
  calibration and translation-invariant moments. Its grouped construction now
  explicitly requires amplitude/noise conditional independence. The focused
  audit, classical CVaR comparators and preferred-view stress results are retained.
- Candidate-only score design uses the three actually fitted Gaussian maps,
  rather than reference-derived alternatives. This still assumes a nominated
  deletion family and tests the complete candidate; it does not certify local
  occupancy or experimentally calibrated density coverage.
- A classical constrained Fisher baseline with independent training/calibration/
  test samples improves a large 10028 perturbation but leaves smaller errors
  and 10076 weak. All 9,450 pointwise difference intervals are retained.
- Equal-noise-budget replica allocation reduces estimated envelopes but often
  enlarges confidence bounds. It does not consistently rescue power and is not
  presented as optimized allocation or repeated-calibration validation.

New simulations are complete. Independent numerical replays and 25 targeted
tests pass, but they do not establish scientific utility or experimental model
adequacy. All noise, transfer and viewing-law bounds are specified assumptions.
The whole-candidate versus local-density distinction and the absence of actual
independent calibration repetitions remain explicit.

The complete published/source evidence is supplied by study identity, including
all unfavorable outcomes and the earlier reviews. The reviewer has read-only
access to exact source/report/JSON/CSV copies and all rendered pages; large NPZ
arrays and third-party publications are not part of this text-tool review.
They must not be described as independently executed or read merely because
manifest hashes are available. Please assess whether any new branch warrants
an ICML contribution, and prioritize essential changes over more similar
simulations. No acceptance conclusion is requested or presumed.
