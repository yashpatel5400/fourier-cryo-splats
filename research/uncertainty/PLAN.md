# Uncertainty-focused research program

Started 2026-09-29. This supersedes the focus of release v0.1.0; that release
remains the record of the earlier reconstruction feasibility study. No new
experimental result or acceptance-level claim is established by this plan.

## User objective and completion criteria

Develop a substantive uncertainty/validation contribution for cryo-EM, informed
by a comprehensive current literature survey; implement and execute demanding
comparisons and failure controls; provide correct supporting theory; completely
rewrite the ICML manuscript; publish reproducible results. Obtain independent
reviews from the actual `claude-fable-5-1` model after a complete submission
candidate exists, address substantive critiques, and repeat. Preserve all
reviews and responses. Never ask the reviewer to manufacture acceptance or
claim that a favorable model review guarantees conference acceptance.

The authenticated Claude Code installation successfully invoked the exact
requested model in an availability-only probe. The returned provider metadata
identifies `claude-fable-5-1`. The first scientific review is now complete; see the current status below.
`caffeinate -im` is running during active work (exec session 14198).

## Work sequence

1. Survey and evidence ledger: uncertainty in density, pose, heterogeneity,
   atomic models and learned priors; validation/FSC and acquisition geometry;
   benchmarks and open questions; general inverse-problem inference.
2. Derive a method with a precise estimand, assumptions, nuisance treatment,
   calibration protocol, and a defensible novelty statement against nearest work.
3. Freeze experiment designs and splits before final evaluation. Pilot studies
   must remain labeled as development; retain failures and method revisions.
4. Numerical theorem checks and synthetic coverage studies under correct and
   deliberately violated assumptions; multiple seeds and Monte Carlo errors.
5. Real-geometry semisynthetic tests with known truth, CTFs, pose perturbations,
   missing views, heterogeneous states, correlated noise and representation error.
6. Real-particle experiments across at least the existing three accessions, with
   source-group splits and held-out predictive checks. Never interpret a real-data
   predictive guarantee as coverage of the unobserved true density.
7. Relevant uncertainty baselines, calibrated and uncalibrated comparisons,
   ablations, efficiency, coverage/width tradeoffs, and negative controls.
8. Write an ICML main paper plus full proofs, survey, and experimental appendix.
   Render and visually inspect the PDF, validate artifacts and reproducibility.
9. Independent Fable 5.1 review and substantive revision cycles; record exact
   model, paper/code hashes, unmodified verdict, and response to every concern.
10. Publish a new version without overwriting the v0.1.0 release history.

## Initial hypotheses, subject to the survey

The central candidate question is whether an uncertainty method can distinguish
random reconstruction variability from bias and lack of identifiable information
under anisotropic acquisition and uncertain poses. Density confidence, observable
prediction confidence, population uncertainty, and map enhancement confidence
must remain distinct. Gaussian Fourier kernels are a computational representation,
not automatically the scientific contribution.

Classical optimal-recovery and bias-aware inference already provide general
linear-functional confidence tools. Any use of these must explicitly credit
that literature; standard Gaussian posterior or conformal formulas are not new
theorems. A cryo-EM-specific contribution must go beyond rebranding them.

## Current status (updated 29 September 2026)

The rewritten ICML manuscript now presents a conditional uncertainty audit in
continuous L2 density space, with explicit nonlinear pose bounds, an integrated
moment refinement, independent numerical checks, and a critical survey of 93
curated candidates. The ledger distinguishes discovery, retrieval and targeted
reading. Forty-eight tests pass. General bias-aware inference, norm duality and
moment inequalities are credited as prior art; scientific novelty and practical
value still require independent assessment.

Completed development covers all three experimental acquisition geometries,
matched covariance/variational/bootstrap controls, grid and continuous audits,
assumption violations, nonlinear feasible adversaries, real-particle Gaussian,
voxel and stock cryoDRGN fixed-pose neural reconstructions, and all frozen
additional-exposure prediction comparisons. Neural prediction wins all six
prespecified paired contrasts on the 4096-particle-per-stack fresh cohort.
This is conditional prediction, not unknown experimental-density coverage.

The original continuous uncertainty validation was frozen at 76dbd53 before
outcomes (96 settings). Its follow-up integrated-moment validation was frozen at
8d66785 after development but before its own outcomes (48 settings, disjoint
particle subsets and new signal/noise seeds). Both are complete with all locked
dependencies unchanged. V1 has minimum in-class coverage 0.993772 over 1536
records; v2 has coverage rounded to one over 768 in-class records. All 72 unique
fits converge. In v2, no-data fallbacks fall from 48 to 16, but correct-sign power
for the in-class cell generators is at most 2.3e-5. Preserve these limitations
alongside the improvements.

The broader-angle development sweep exposed no-data fallback in every original
1/2/5-degree setting. The refinement reduces fallbacks from 18 to 10 across 30
settings. All 108 optimized continuous feasible adversaries stay below its
upper bound, but the best attained biases are only 0.16--0.57 of that bound.
Conservatism remains. Higher-frequency continuous solves are complete for all three geometries;
runtimes reach 1223 seconds per target and reference sign power remains poor.
A classical fixed-length lower bound is within a factor 1.20 of every reported
fixed-pose interval; it does not establish the larger pose bound's sharpness.
Real-image background diagnostics expose unresolved dependence after scaling.

The first authentic Fable 5.1 review is complete: reject, confidence 4/5, not a
strong ICML contender. Its original prompt, source/PDF hashes, provider output
and verdict are preserved under reviews/round-01. The reviewer found the checked
mathematics correct but practical calibration, usefulness and pose scaling
insufficient. The finite revision plan prioritizes continuous pose-aware weights,
one documented experimental calibration attempt, a matrix-free pose scale test,
and realistic translation/CTF sensitivity. These are substantive open tasks.
The next major decisions should follow this assessment, rather than accumulating
more similar oracle-bound simulations. End-to-end pose/noise/class calibration,
biologically useful feature resolution and the strongest matched external
baselines remain material scope questions. The goal is active, not complete.

The first revision now includes a tested matrix-free pose operator, a numerical
failure-probability budget, a small independently checked pose-aware optimizer,
and completed experimental-noise/CTF sensitivity attempts. The latter are
negative usefulness results, not established experimental calibration. The
10,000-particle scale run and longer three-geometry pose-aware fits are active.
Their final width, adversarial sharpness and reference-power results should
guide the next methodological decision. The signed density class remains too
broad for the desired fine-feature interpretation; any support/positivity or
energy refinement must be explicit, physically motivated and independently
validated rather than chosen to turn a negative outcome positive.
