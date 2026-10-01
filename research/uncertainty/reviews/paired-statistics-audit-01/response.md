# Response to the focused paired-statistics audit

The exact requested `claude-fable-5-1` returned one preserved text report.
It checked the formulas favorably but did not issue an acceptance verdict.
All supplied source hashes were unchanged during the invocation. It did not
execute code or inspect numerical arrays. The complete report remains unchanged.

- **P1/P2:** implemented the spectral two-sided growth diagnostic, kernel
  nonnegative covariance projection and all 440 real coordinates. The earlier
  compressed results remain. An SVD checks all 18 feature matrices: each has
  full column rank at relative threshold 1e-8. Rank-32 discrimination retention
  is .819/.973/.828, rather than a claim that preserving signal energy preserves
  every detectable alternative. The full method still has optimization gaps.
- **P3:** checked 10,000 new orientations per stack, then refit on all 14,160
  views and checked a third 10,000-view set. Finally 80 continuous rotation/shift
  searches on the surviving stack find 42 positive exact-forward violations.
  Four worst cases are independently verified by physical-cell exponential
  sums within 8e-14. The saved directions cannot be promoted to continuous-pose
  tests. No sampled maximum is called a continuous certificate.
- **P4:** all twelve common-shift cases finish at SD 0,.5,1,2 downsampled
  pixels. Alternative moments integrate the Gaussian law exactly using both
  complex second moments, checked by independent Gauss-Hermite quadrature.
  Shift sensitivity and every fresh-view failure are retained. These checks
  do not solve arbitrary continuous translations or relative frame motion.
- **P5:** the new report distinguishes the diagonal near-membership result
  on 10076 from open discrimination on 10049. It does not use the 128-particle
  expected score as a power cutoff. Stored normal approximations at larger
  sample sizes are labeled as approximations and were superseded by actual
  pose counterexamples; no measured rejection-power claim is made.
- **P6:** min/max mixture-mass LPs on the same 4,160-view catalog give
  [.964499,1.154280] for half removal and [1.043575,1.146993] for full removal.
  Half removal can match at unit mass by interpolation; full removal needs
  amplitude at least about 1.0216 in this finite catalog. The earlier arbitrary
  LP vertex is not a canonical amplitude and these are not continuous extrema.
- **P7:** added primal/dual consistency and raw dual-feasibility assertions,
  plus a zero direction for near-membership. Historical weights and records
  are preserved. The full-frequency bounds remain valid for any nonnegative
  mixture without claiming that solver convergence proves global test optimality.
- **P8:** all three frequency sets have 220 unique nonzero frequencies, no
  opposite pairs and no Nyquist coordinates. Identity covariance is stipulated
  in these oracle screens, not inferred from experimental spectra. Large null
  simulations and paired-empty-ice calibration are deferred because the key
  continuous-pose condition has concrete counterexamples.
- **P9:** a new scalar evaluator adds the concave tangent endpoint bound to
  cover an inexact stationary root, plus an explicit diagnostic rounding pad.
  All 36 enlarged-cone bounds are recomputed. The 10076 full-removal unit-noise
  upper for 128 particles is now 6.411e-10, not the old unpadded 7.55e-28.
  Frozen outputs stay public; ordinary floating point remains unvalidated.
- **P10:** added a proper Gaussian moment-domain eigenvalue check, singular
  covariance and near-boundary tests, an off-frequency translation check,
  an independently optimized spectral scalar reference and a direct Gaussian
  log-variance simulation. The saved-array replay is explicitly independent
  of the LP, not an independent dual implementation. A local directional
  search is never labeled its global maximum. Unequal-mean physical exposure
  effects remain a limitation rather than an asserted calibrated case.

The [complete outcome report](../../paired-power-v1/METHOD-GATES-RESULTS.md)
retains all cells and numerical corrections. These steps improve diagnosis
and falsify premature validity claims; they do not complete a useful, calibrated
new uncertainty method. All three full ICML reviews still reject.
