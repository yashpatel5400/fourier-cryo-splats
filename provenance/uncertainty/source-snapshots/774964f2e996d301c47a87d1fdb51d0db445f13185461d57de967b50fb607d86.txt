# End-to-end local-pose interval study, version 1

1 October 2026 UTC. Motivated by full review 2 R9/R10/R11. Sources and this protocol must be committed before any study interval/coverage output. The preceding local-alignment calibration is separate, with seed root 261001 and 128 replicate datasets per geometry. Early calibration pose errors have been inspected; no study coverage exists at declaration. No parameter may be chosen from the following 200 test replicates.

## Population and replication

Use all three existing radius-five acquisition geometries (128 particles, seed 609311), the saved pilot-registered unit-L2 reference generators and independent pilots. Test exactly 200 independent replicate datasets per geometry, seed root 261002. Truth and geometry are fixed within a geometry; new Gaussian noise and local-pose initialization are drawn on every replicate. The 10076 Class A generator is a homogeneous simulation, not truth for the heterogeneous experimental stack. These are three generators/acquisition designs, not 600 independent biological systems.

Generate alignment image Y_A and independent image Y_B at the same noiseless signal and prescribed scalar noise. For each replicate fit poses to Y_A using both the oracle-reference and independent-pilot templates, with exactly the frozen calibration optimizer/initialization distribution. These controls share images and initializations. Do not retry or select a template. Dephase Y_A and Y_B by their estimated detector shifts and construct a new continuous Fourier Gram using their estimated rotations. Recompute all estimator weights. This is a local-alignment pipeline; density/template learning and global orientation discovery are not repeated.

Two fixed broad targets: a central Gaussian average and a signed contrast at z=+0.08 and −0.08 of the field, both SD 0.07 field. This first study addresses R11. It does NOT satisfy the separate radius-12 / ≤20-A usefulness criterion R2. Each procedure/target reports marginal 95% intervals; no simultaneous 95% claim over the whole table or minimum across methods is made.

## Calibration rules, fixed before study coverage

For each geometry/template, use all 128 calibration replicate pose arrays. Define the score for a whole dataset as the maximum over particles of sqrt((rotation error in degrees)^2 + (translation error in Angstrom / 0.5)^2). Select the largest of the 128 calibration dataset scores. Under exchangeability of calibration/test datasets with this SAME fixed truth/geometry/template/optimizer, the chance a fresh score exceeds that maximum is at most 1/129. This is a simulation-only known-truth tolerance construction. It is not experimental pose calibration and is not a guarantee uniform over a new population of molecular densities.

The robust audit uses joint radii angle=score degrees and shift=0.5*score Angstrom, with measurement alpha=.04. For the independent Y_B control, conditional fixed-design bias coverage plus the marginal tolerance event gives miscoverage at most .04+1/129<.05. It does not justify same-image Y_A inference. Save realized score/bound failures and raw widths even when a no-data fallback is used.

The mixed-pose alternatives use these same bounded radii, independent normalized radius 1, and common-mode radii 0 and 0.1. They assume independent centered residual errors conditional on weights/design. Estimated poses do not establish that assumption, even with independent Y_B. Their plug-in coverage is an empirical test, not a consequence of the conditional theorem. The 0.1 option is retained separately; it is not selected per outcome.

For each local-Gaussian-pose baseline, set the five coordinate SDs to the calibration root-mean-square error per coordinate over all particles/replicates, without subtracting the empirical mean. Save means as a bias diagnostic. This is a declared linear Gaussian working model, not calibrated nonlinear nuisance uncertainty.

## Nine procedures, both observation controls

1. Fixed-pose folded-normal audit.
2. Fixed-pose sum-width audit (bias plus normal quantile times SD).
3. Deterministic first-order-plus-second-remainder pose audit, using the simulation-calibrated joint radius and alpha=.04.
4–5. Mixed centered independent/common-mode intervals, common radii 0 and 0.1.
6–9. Continuous generalized Gaussian-prior intervals at tau=B=2 and tau=B/2=1, both fixed pose and first-order Gaussian pose marginalization.

The first five share the fixed-pose sum-objective weights; their raw uncertainty bounds are distinct. Choose a no-data fallback at the pilot centre only if that raw bound exceeds B||ell||. Save every raw and selected interval, and whether the fallback was used. Bayesian working-model intervals are reported raw. Radius-five quadrature order 40, preconditioner rank 1024 (selected by a no-pixels/no-coverage numerical preflight), fixed-audit maxiter 100 and relative gap .005; Gaussian CG tolerance 1e-10, maxiter 2000. All solver diagnostics and failures are retained. These are ordinary floating-point calculations with analytic integration control, not validated roundoff enclosures.

## Reporting and failures

Save a JSON outcome and weight/pose/observation NPZ for every replicate, both templates, both targets, both image controls, every method. Preserve nonconvergence and exceptions; a failed replicate may not be omitted from the planned denominator of 200. Report empirical raw/selected coverage, exact binomial 95% Monte Carlo intervals, widths in feature units and relative to no-data, fraction with relative width<0.5, sign-exclusion probability, fallback rates, observed pose errors, dependence/mean-error diagnostics, solve gaps and runtime/memory. Images/templates/methods/targets within a replicate are paired, not independent replications. No optional stopping. Summaries cannot be marked complete until every prescribed replicate exists.

The independent-image robust theorem is narrow and simulation-assisted. A favorable mixed-model result would still not establish its centering/independence premises, experimental covariance transfer, an appropriate physical density class, or publishable novelty. Unfavorable coverage or vacuity remains a central outcome.
