# Response plan to round 1 (not a completed response)

The authentic review rejects the candidate. We accept its central assessment:
conditional correctness and extensive implementation checks do not establish
practical usefulness. No favorable review or scientific completion is claimed.
The original review and exact candidate packet remain unchanged.

## Essential revision work

- **R3:** Optimize continuous estimator weights with pose uncertainty included.
  First complete a bounded bandwidth-choice development probe on existing
  geometry/noise, then develop a shared-field convex pose objective. A simple
  sum of independent frequency bounds is much looser; it is not an adequate
  replacement. Decide gauge explicitly: current targets are in the declared
  external coordinate frame, so a coherent rotation belongs in the stated
  sensitivity class. A gauge-quotiented alternative requires changing both
  target and nuisance class, not deleting unfavorable adversaries.
- **R4:** Develop a matrix-free field/adjoint operator and a defensible spectral
  upper bound. A Lanczos Ritz value alone is a lower bound, so it must not be
  substituted for the current spectral upper audit. Demonstrate measured time
  and memory at 10,000 particles, retaining any loss of useful precision.
- **R1/R8:** Attempt one experimental uncertainty analysis on 10028 using
  held-out exposure noise diagnostics, documented pose information, declared
  density assumptions, realistic translations and CTF sensitivity. Estimated
  nuisance bounds must not be called calibrated without justification.
  A deposited-map comparison is an approximate-reference diagnostic, not a
  density-coverage label; retain any vacuous result.
- **R2/R7:** Test whether the resulting method is useful at a higher band and
  structural scale, with a defensible support/positivity class if applicable.
  Neither a relabeled classical theorem nor more similar simulated cases
  addresses the novelty/usefulness objection.
- **R5/R6:** Reframe around continuous Fourier-slice auditing and its pose
  component; make simulated-noise/known-bound conditions explicit in the
  abstract. Lead with widths and power. Move reconstruction prediction to
  supplementary context. Retain all original frozen outcomes and Monte Carlo
  checks as historical verification rather than deleting inconvenient records.

## Corrections, clarifications and smaller items

- **M1:** This is a public non-submission preprint, not an anonymized ICML
  submission. Replace the misleading anonymous placeholder; prepare an
  anonymous version only if submission is later requested.
- **M2:** Historical survey counts are now explicitly marked as historical.
  The current curated count is 93; this is not a full-reading count.
- **M3:** Verify the raw bootstrap count rather than choose a rounding to match
  either prose or the review.
- **M4:** Remove the asymptotic critical-value shortcut in a new revision
  implementation while preserving all frozen source dependencies and their
  original outputs. Audit old numerical consequences separately.
- **M5/M6/M7/M11:** State endpoint checkpoint selection, physical shift units,
  what each continuous adversary represents, and the supplied simulation
  assumptions directly.
- **M8:** Check the primary Low/Cai--Low papers and cite their actual relevant
  results; do not attach a citation based only on a review suggestion.
- **M9:** Keep the user's requested comprehensive survey as a standalone
  research artifact; shorten the submission-style appendix as appropriate.
- **M10:** The original evaluators abort on undercoverage before preserving the
  failing case. No such abort occurred in either complete study, but the
  failure-recording policy is flawed. Future evaluators must persist the case
  and diagnostics before flagging failure. Frozen evaluators will remain
  reproducible at their original version, with this flaw disclosed.
- **M12:** The reviewer saw source, not pages. The rendered round-1 main text
  ends on page 7, before references; it does not exceed eight pages. This
  factual correction does not remove the legitimate focus/clarity criticism.

Other factual scope clarifications: spatial annulus correlations do not prove
that selected Fourier noise coordinates violate the assumed covariance, since
the background may contain signal and preprocessing effects. The paper
explicitly refrains from that inference. Nonetheless the absence of validated
experimental whitening is real and remains unresolved. The higher-band and
1,024-particle studies exist, but neither supplies the missing large-scale
pose-aware usefulness result.

No item is marked resolved merely because this plan names it. A subsequent
review will receive the unaltered prior review and an evidence-linked response.
