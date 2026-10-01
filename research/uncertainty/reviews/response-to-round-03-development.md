# Response to full independent review 3 — active development

1 October 2026 UTC. The unchanged review recommends **reject**, confidence
4/5, and explicitly says the work is not a strong ICML contender. All 26
rendered pages were visible. The exact provider model was `claude-fable-5-1`;
all 2,437 evidence copies were unchanged. The reviewer read code and records
but executed no tests, opened no numerical arrays and read no third-party
papers. Favorable mathematical checks do not resolve the rejection.

## Accepted priorities

- **R13:** reanalyse every saved trial for standardized alignment bias, RMSE
  against the pilot, noise-only and reduced-radius coverage, paired image
  differences, and an added true-pose arm. The [post hoc protocol](../REFITTING-REANALYSIS-PROTOCOL.md)
  distinguishes new analyses from the original frozen study.
- **R14:** the nonlinear remainder is not the sole obstacle. Even the
  first-order tail alone exceeds no-data width. Compute exact class-worst-case
  bias at realized poses before proposing a new construction; do not launch
  another series of Taylor-bound refinements or unchanged trials.
- **R11/R9:** the truncated coarse maximum-likelihood optimizer fits noise.
  A subsequent simulation requires an estimator that improves initialization
  in a pre-study gate, true-pose/noise-only controls, and a boundary generator.
  The conditional centering premise remains unestablished.
- **R7/R15:** make the ridge-path equivalence explicit, investigate folded-width
  optimization, and disclose that the class radius has no experimental
  calibration. Existing operator integrals alone do not establish novelty.
- **Reporting:** add class-worst-case Gaussian coverage at both prior scales,
  use the frozen prediction outcomes, correct the bibliography mismatch after
  checking the primary sources, and replace redundant coverage tables with
  the alignment-bias figure while retaining every row in the repository.
- **R1/R2/R8:** experimental nuisance calibration and practical fine-scale
  inference remain open. The new paired-power candidate is unreviewed and has
  only primitive mathematical tests, not an established substitute method.

## Two factual distinctions

The review says the simulation draws true poses around supplied initializations.
The actual runner fixes the true poses (`signal` is fixed), then draws a noisy
initial rotation and shift. A latent-pose posterior model is a possible new
study, but is not the probability model of these archived trials. This
correction does not rescue the centering assumption or the optimizer.

Seven of 200 scores falling outside one calibrated ball gives an empirical
3.5% conditional failure estimate, not a known conditional population failure
probability. The 1/129 tolerance guarantee is marginal over calibration batches;
neither this distinction nor the observed inclusion rescues the vacuous bounds.

## Status

The diagnostic analysis is being implemented. No new full review or acceptance
assessment has occurred. All three full reviews remain unchanged and public.
