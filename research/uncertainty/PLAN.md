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
identifies `claude-fable-5-1`. No scientific review has yet been requested.
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

## Current status

Literature audit underway. Main methodology is not frozen. Small physical
development experiments and numerical/conic checks have run; their limitations
and a necessary method revision are recorded in DEVELOPMENT-LOG.md. Three EMDB
reference volumes were downloaded as semisynthetic signal generators. Whole
micrograph/film splits are now available for all three original datasets.
For EMPIAR-10076, the deposited Frealign FILM field resolves the misleading
sequential micrograph labels in its STAR file. The main manuscript has now been
completely rewritten as an explicitly labeled uncertainty development draft.
The first ambient-space audit and full-grid solves cover all three acquisition
geometries; their known-map coverage is conditional on fixed poses, Gaussian
noise, a finite supported grid and a prescribed norm class. Final validation,
method integration and independent scientific review remain outstanding.

The next methodological priority is combining the ambient density-space audit
with useful nuisance bounds, rather than reporting finite-dictionary pose
coverage and full-grid fixed-pose coverage as though they were the same result.
Fresh confirmatory data are still required after development freezes.
