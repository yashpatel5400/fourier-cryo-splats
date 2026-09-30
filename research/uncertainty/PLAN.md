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

## Current status (30 September 2026)

The full round-1 Fable rejection remains the latest acceptance assessment.
The survey contains 106 curated candidates, with targeted reading distinguished
from retrieval. The 10,000-particle audit, expanded Fourier variational baseline,
three rounds of alternating pose/density ambiguity, and enclosing-domain
remainder refinements are complete. Their improved numerical bounds have not
established fine-scale experimental usefulness. Numerical tests verify
implementation; their count is not evidence of scientific acceptance.

The eighteen-case pose-aware comparison is complete. Its two-degree results
do not meet the reviewer's usefulness or bound-tightness criteria. All twelve
locked pilot-selected fixed-pose fits are complete; their nonlinear audits and
experimental applications are running. The latter uses three dataset processes, each with two
FINUFFT threads, to exploit local resources without changing the protocol.
A fresh noise-calibration cohort reserves 128 previously unused exposures per
stack. All twelve estimators and 116 files were frozen and published in commit
60efd9b before any reserved pixels were accessed. All three stacks are now
downloaded and verified. The coordinator awaits the complete original
twelve-feature application before running the frozen recalibration. The old inference
images are still development data, and pose/density assumptions remain open.

The metadata availability audit finds zero-valued pose/shift ESS fields in all
three archived CS files; these are not zero-error bounds. In particular, the
single-density inference model is not experimentally validated for heterogeneous
EMPIAR-10076. R1, R2 and R7 remain material scientific objections. The goal is
active, with no claim of an acceptance-ready submission.

The matched pilot-target Fourier baseline now retains both its original prior
sweep and a declared post-outcome broader-prior sensitivity. Broader priors
produce favorable reference results and are not omitted. A same-weight continuous
sign/total-norm ablation is also complete: stronger sign constraints exclude
the processed reference, and all one-degree cases still lack sign power.
These are development comparisons, not new frozen density calibration.

A completed fixed-weight cubic audit reduces one selected one-degree width
from 0.473 to 0.252 of no data, without useful sign power. Its separately
declared cubic weight fit completes in 7,191 seconds with original nominal-Gram
coordinates. Its relative width is 0.180 and minimum reference sign power 0.00654;
the iteration limit and 0.940 relative surrogate gap preclude a convergence claim.
Two-tolerance final-weight checks are stable, without validating all floating-point
error. Focused Fable audits support the real-arithmetic argument but do not
constitute an acceptance review.
A separate pose-penalty coordinate metric passed its compute-only probe and a
small independent conic optimization check. The original gap triggers the
previously declared empirical follow-up using that metric; it is now running
from the same original weights, with no asserted convergence improvement.

The local Gaussian pose-marginalization baseline is complete: 48 solves across
the same twelve targets and two broader priors, with 384 analytic conditional
reference checks. It preserves all fixed-pose controls and both favorable and
unfavorable outcomes. Its small width inflation is specific to the pilot
linearization, whose nonlinear error and omitted density/pose products remain
limitations. This supplies an additional matched control, not new experimental
calibration or a full nonlinear pose posterior.

A native RELION 5.0.1 CPU probe completed in 14.794 seconds on 256 old pilot
particles. A separately published three-stack unknown-pose reconstruction batch
is running: pilot-only VDAM initialization, then refinement on the same exposure
halves as the existing neural comparison. It has explicit per-stage wall limits,
retains failures, and has no converged baseline result yet. Global map alignment
is selected only against the old pilot, with all 48 local starts retained; two
analytic-phantom checks pass. This is preparation for an external reconstruction
comparison, not a new density-uncertainty guarantee.

## Historical status (29 September 2026)

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

### Subsequent revision checkpoint

The 10,000-particle matrix-free scale test is complete (1.95 GB peak memory).
The curated survey now contains 98 candidates, not 98 full readings. Four
focused mathematical audits and their responses follow the unchanged full
round-1 rejection. Pilot-specific analytic moments reduce a two-degree
full-weight example to 0.510 of no-data width, without useful reference sign
power. A 10 Å Gaussian-width, 1,024-particle fixed-pose fit is precise, but its
one-degree pose audit returns no data; the cubic remainder dominates. The
three-stack pose-aware grid is not complete.

The newest prospective component constructs two-pose/density pairs to lower-
bound the width any deterministic-length interval needs on the declared
class. It is a classical testing specialization intended to diagnose the gap
between intrinsic ambiguity and conservative upper bounds. Independent conic
checks pass and the initial 30-case three-stack grid is complete. A dedicated
Fable mathematics audit confirms the argument and identifies numerical and
recording gaps. These are fixed; the unchanged grid reruns with additional
witness records and error guards; all 30 revised cases are now complete, with
fixed-pair brackets below 0.5 percent. One hundred and two tests pass; that does not resolve
the scientific usefulness or experimental-assumption objections. No full
second-round acceptance assessment has occurred, and the goal remains active.

The current enclosing-domain remainder refinement retains Fourier cancellation
within each particle. The ball and expanded-cube probes reduce the fine-target
cubic term from 56.493 to 32.623 and 21.110, respectively; the latter interval
has width 0.751 of no data but still no useful sign power. Independent integration
and nonlinear tests pass (109 tests total). This moves the numerical bound but
does not settle R1/R2/R7 of the full review. The next meaningful objective remains
useful precision with defensible experimental nuisance/noise/class assumptions,
rather than the number of tests, probes or appendix pages.
