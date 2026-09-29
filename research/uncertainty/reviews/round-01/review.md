**Verdict up front: reject.** Confidence 4 of 5. The current manuscript is not a strong contender for ICML acceptance. The mathematics I could check by hand is correct, the provenance and freezing discipline is unusually good, and unfavorable results are reported plainly. But the delivered method is a specialization of classical bias-aware inference whose intervals are informative only in a toy regime, whose required inputs cannot yet be obtained from experimental data, and whose one genuinely new component, the continuous nonlinear-pose bound, is loose by a factor of two to six exactly where practical pose errors live.

**What I could not verify.** No tools were available. I did not execute code or tests, did not see any figure, and did not inspect per-case JSON files, lock manifests, or git history beyond the supplied summaries. Development tables in the appendix, including the ambient-audit and ambient-pose tables, have no supporting result files in the packet. I read no third-party papers and rely on the project's own annotations plus my background knowledge of the cited work. All formula checks below are by hand.

## 1. Contribution and closest prior art

**What the paper actually delivers.** For a fixed linear density functional, an affine estimator centered on an independent pilot, and a declared class of densities and nuisance parameters, the paper computes a bias-aware confidence interval whose worst-case bias term is evaluated in a class larger than the fitting representation. The concrete new ingredients are these.

- **Analytic continuous audit.** The squared L2 norm on the unit cube of the realified Fourier-slice adjoint field is a sinc-kernel Gram, and the Gaussian target pairings are truncated-normal integrals with complex arguments. This gives the fixed-pose bias bound in continuous L2 without a voxel grid or bandlimit.
- **Matrix-free quadrature Gram** with an explicit Gauss remainder that keeps both primal and dual bounds valid in real arithmetic.
- **Bounded nonlinear pose corrections** for the Fourier-slice operator: first and second pose-derivative fields from twenty Fourier moments per particle, a shared-density spectral bound, a uniform cubic remainder, and an integrated-moment refinement of that remainder.
- **Frozen simulation studies** on acquisition geometries from three EMPIAR stacks, plus assorted development controls.

**Closest supplied prior art.** Armstrong and Kolesár supply the bias-aware critical value and the uniform-over-class coverage argument. Donoho supplies optimal recovery, the dual characterization, and the fixed-length efficiency factor near 1.2. Ullrich et al. supply differentiable Fourier-slice reconstruction with explicit density uncertainty and the variance-versus-bias distinction. Rangan et al. supply pose–density soft modes. Lai et al. supply scalar confidence sets in a cryo-EM setting. Batlle et al. supply simultaneous constrained inverse-problem regions. The manuscript credits all of these. Its own description, "explicit continuous Fourier-slice specialization," is accurate. Nothing in the statistical construction is new. The novelty is in the cryo-EM-specific integrals and pose constants and in the implementation.

**What the paper does not deliver.** It does not deliver an uncertainty method that can be run on experimental particles with a guarantee, because the noise law, the pose radii, and the density radius are all supplied by assumption and the paper's own diagnostics show the experimental noise assumption is violated. It does not deliver informative intervals for features finer than roughly 15 Å, or at rotation budgets of one to two degrees on two of three stacks.

## 2. Technical checks

I re-derived each claimed identity and bound. Results are summarized here; details follow for the ones that needed care.

| Item | Status | Note |
|---|---|---|
| Theorem 1 conditional coverage | Correct | Standard folded-normal argument. Requires weights and fallback decisions independent of inference noise, which the frozen pipeline satisfies. |
| Adjoint field and sinc Gram | Correct | Realified pairs give the field as the real part of a complex sum; the cross terms carry the sign-flipped kernel; the factor one half is right and matches the code. |
| Truncated Gaussian transform | Correct | Completing the square gives the complex-argument normal CDF form. The squared target norm uses the product-of-Gaussians identity with variance halved. |
| Continuous dual and feasible lower bound | Correct | Weak duality holds. Padding the residual norm upward and the observation norm upward keeps the dual point feasible. |
| Gauss–Legendre constant and E_r | Correct | I reproduced the order-40, radius-5 remainder of about 3e-23 by Stirling's formula. The realification factor of root two is justified. |
| Polynomial-exponential quadrature error | Correct | Leibniz bound, product frequencies at most twice the maximum coordinate, monomial factors bounded by one on the cube. |
| Grouped spectral bound and Frobenius pads | Correct | Elementary operator-norm inequality; integration pads enter conservatively. |
| Pose derivative fields | Correct | Row-vector rotation convention is consistent between derivatives, generators, and tests. The rotation Hessian of the phase is correct, and the scaled symmetric vectorization preserves the unit-ball norm. |
| Cubic remainder | Correct | Third derivative of the phase exponential has the three stated terms. |
| Integrated moment constants | Correct | Uniform even moments 1/12, 1/80, 1/448 are dominated by matched-variance Gaussian moments; Minkowski and Hölder give the stated envelope. |
| Fixed-length lower bound and 1.198 factor | Correct but tautological | This is Donoho's fixed-length affine efficiency bound. The empirical ratio only restates that the optimizer converged. |
| Gaussian dictionary Gram | Correct | Even and odd parts are orthogonal, giving the two blocks. |
| Chi-square noise-scale bound | Correct, textbook | Independence of calibration coordinates is essential and, per the paper's own control, fails under equicorrelation. |

**Two subtleties worth stating explicitly.** First, the pose Taylor remainder is charged per particle by the triangle inequality while the polynomial part is bounded jointly. That is valid. Second, the coherent rotation used as the boundary adversary applies one identical rotation to every particle. That is a global gauge rotation of the reconstruction frame, which the paper elsewhere says must be fixed or included in the target class. The worst-case pose bound therefore contains a mode that is arguably a coordinate choice rather than an alignment error. See R3.

**Consistency of quoted numbers.** Every headline number I could check against the supplied summaries matches, including minimum in-class coverage of the two studies, all width ranges, the no-data fallback counts, the sign-power maxima, the remainder reduction range, the adversary ratios, the higher-band timings and widths, and the background pair-product range. Three numbers could not be checked: the development pose audit's feasible-bias ratio range, the enrichment convergence counts, and the ambient-audit table.

## 3. Experimental audit, concerns, revision plan, and verdict

**Frozen versus development separation.** The separation is real and documented: protocols, seeds, lock hashes, and disjoint particle positions are all recorded, and the moment refinement was explicitly developed after seeing development failures and then re-run on disjoint particles. Two caveats. The frozen studies reuse the same three stacks, pilots, targets, noise convention, and class as development, so nothing new about generalization across acquisitions was tested; replicate widths vary by under five percent and add little. And the evaluator raises an assertion and stops on any in-class undercoverage, so a failure would have produced an incomplete study rather than a recorded failure, contrary to the protocol's "retain every case."

**What the frozen coverage numbers mean.** Coverage in every record is computed analytically from a deterministic bias and a known Gaussian noise scale. The ten thousand draws per case only check arithmetic. Given Theorem 1, in-class coverage at or above the nominal level is guaranteed whenever the implementation is correct, so the studies are implementation validation. The informative outputs are widths and sign power, which are summarized below.

| Setting | Feature scale | Relative half-width | Sign power on map-like signals |
|---|---|---|---|
| v1, 0.1°, broad | 16–34 Å | 0.11–0.17 | 1.0 for central averages; near zero for contrasts |
| v1, 0.5°, broad | 16–34 Å | 0.27–0.46 | mostly near zero |
| v1, 0.1°, fine | 7–15 Å | 0.76–0.78 | near zero |
| v1, 0.5°, fine | 7–15 Å | 0.89–1.00 | zero |
| v2, 1°, broad | 16–34 Å | 0.30–0.50 | at most 2.3e-5 |
| v2, 2°, broad | 16–34 Å | no data on two stacks; 0.87–0.99 on the third | zero |

**Regime and assumptions actually tested.** All uncertainty experiments use 128 particles, Fourier radius five, simulated white Gaussian noise with a known scale, exact CTFs, a density radius of two around unit-norm functions, and a translation budget that is negligible in physical units.

| Stack | Field | Radius-5 band limit | Shift budget |
|---|---|---|---|
| 10028 | 482 Å | 96 Å | 0.20 Å |
| 10049 | 236 Å | 47 Å | 0.10 Å |
| 10076 | 419 Å | 84 Å | 0.17 Å |

A radius of two around unit-norm functions is the trivial maximum distance, so the class carries essentially no information from the pilot beyond its center. A shift budget of a tenth of an ångström is far below any realistic alignment uncertainty. The rotation budgets that produce usable widths, 0.1° and 0.5°, are below typical angular accuracy for these specimens and are uniform over particles, which no alignment procedure provides.

**Baseline fairness.** The mathematical controls are honest and correctly matched. The most informative comparison is the exposure-group bootstrap, which on the deposited-map signal attains near-nominal coverage with widths two to four times narrower than the audit. The paper concedes this. No external uncertainty method is run on the same estimand; I do not require one, since the supplied methods answer different questions, but the consequence is that the empirical case for the audit rests on adversarial signals that no practitioner would expect to encounter. The neural-versus-classical prediction study is a fair prediction comparison but is unrelated to the uncertainty claims.

**Claims not supported by the evidence.** The abstract's "144 audit settings and 3,456 signal/scenario records" presents software validation as evidence of calibration. "Within a factor 1.20 of the optimal deterministic length" is a restatement of a classical bound and of solver convergence, not a finding. The introduction's "retain Gaussian computation while auditing feature uncertainty in a larger class" is not what the frozen pipeline does; the continuous solver optimizes directly in the continuous class and never uses a Gaussian dictionary. The claim that the pose bound "does not expose failure at the tested severities" for the fourfold exceedance control shows conservatism, not robustness, and should not appear near the word "covered."

**Major concerns.**

**R1. The intervals cannot be produced for experimental data.** Known white Gaussian noise, exact CTFs, uniform pose radii, a density radius, and an independent design are all assumed. The paper's own background diagnostics show lag-one pixel correlations of 0.2 to 0.3 in all three stacks, its stress tests show correlated noise causes severe undercoverage, and its chi-square calibration control fails half the time under moderate equicorrelation. Consensus poses were refined against all particles. Change needed: an end-to-end run on at least one stack with a noise covariance estimated from held-out exposures, a pose radius from an external or held-out source, a stated physical prior for the density radius, and validation of the resulting interval against the deposited high-resolution map treated as approximate truth with explicit caveats. Acceptance criterion: a nonvacuous interval for at least one named feature that covers the deposited-map value, with all three inputs documented. Status: fatal to any applied-method claim; fixable only by adding a substantial new component. If the outcome is vacuous intervals, that is a legitimate finding, but then the paper is a negative result and should be written as one.

**R2. Usefulness is not demonstrated in any regime that matters.** The band limit of roughly 50 to 100 Å, 128 particles, and a class of signed unit-norm functions with radius two make the feature-level guarantees uninformative for structural biology, and the bootstrap dominates on realistic signals. Change needed: a study at radius twelve or higher, on at least a few thousand particles, with a narrower and physically motivated class such as positivity or a support mask, and reporting of sign power on map-like signals as the primary usefulness metric. Acceptance criterion: an interval that excludes zero for a known feature at 10 Å scale or finer under a rotation budget of at least one degree. Status: fatal to a usefulness claim; potentially fixable, but only jointly with R3 and R4.

**R3. The continuous pose bound is a post-audit of the wrong weights and is loose where it matters.** Feasible adversaries reach only 0.16 to 0.57 of the integrated bound at 0.5° to 2°, and at 2° a feasible bias alone would support widths of 0.16 to 0.27 of no data against upper bounds of 0.86 to 1.46. In the finite-voxel development, optimizing weights with curvature retained gave 0.42 to 0.53 at 2°, and the shared-density post-audit halved that again. The continuous pipeline never optimizes with pose terms. Separately, the coherent-rotation adversary is a gauge mode; the paper should decide whether global rotations belong in the nuisance class and, if not, remove them from both the bound and the adversaries. Change needed: optimize continuous weights under the pose-robust objective, which is already a sum of Euclidean group norms compatible with the existing solver, at 0.5°, 1°, and 2°, and report the feasible-over-upper ratio. Acceptance criterion: ratio at least 0.5 at 2° and widths below 0.5 of no data on all three stacks. Status: fixable and the most valuable single experiment.

**R4. The pose audit does not scale.** It stores twenty fields per particle at all quadrature nodes and forms a dense Gram over all twenty-n columns. Storage is 725 MB at 128 particles and grows quadratically, so ten thousand particles would need terabytes as implemented. The fixed-pose matrix-free solver was run at 1,024 particles but the pose audit was not. Change needed: a matrix-free spectral bound via Lanczos on the field operator with NUFFT matvecs, and a demonstration at ten thousand particles or more with a stated memory budget. Acceptance criterion: audit completes at that scale within the budget with the same guarantees. Status: fixable in principle, undemonstrated.

**R5. The frozen studies are presented as calibration evidence rather than implementation validation.** Coverage is analytic and tautological given the theorem; the correctness criterion cannot fail unless there is a bug; "coverage rounds to one" at one to two degrees reflects conservatism, not calibration. Change needed: reframe both studies as verification plus width and power outcomes, remove record counts from the abstract, drop the redundant Monte Carlo draws, and report widths and sign power per stack as the primary tables. Acceptance criterion: no coverage number is presented as empirical evidence of calibration on unknown density. Status: fixable by rewriting.

**R6. Title, framing, and the prediction study do not match the delivered method.** The Gaussian dictionary and the enrichment scheme appear only in the development appendix, the enrichment converged in eight of twelve cases and is slower than the matrix-free comparator, and the fresh-exposure neural prediction study has no bearing on density uncertainty. Change needed: retitle and reframe around the continuous Fourier-slice audit and pose bounds; move the prediction study to supplementary material or drop it; keep the dictionary-exclusion demonstration, which is the paper's clearest motivating result. Status: fixable by rewriting.

**R7. Novelty is incremental relative to the supplied prior art.** The statistical machinery, the duality, the efficiency factor, and the noise-scale bound are all classical, as the authors say. What remains new is a set of integrals and Taylor constants for one operator. That can justify a paper only if R2 and R3 show the specialization produces something usable. Status: not fixable by wording alone; contingent on R2 and R3.

**R8. The nuisance model omits realistic components.** Translations are effectively fixed, pose radii are uniform across particles, CTF parameters are exact, and per-particle scale and B-factor errors are absent. The low-bandwidth defocus control passes only because CTF barely matters at 50 to 100 Å. Change needed: per-particle pose radii and a physically plausible translation budget of at least half an ångström in the main study, plus a CTF-error term at the higher band. Acceptance criterion: widths reported under those settings. Status: fixable; likely to widen intervals further, which must be reported.

**Minor concerns.**

- **M1.** The manuscript lists "Anonymous Authors" but includes a repository URL containing the author's handle. This breaks double-blind review.
- **M2.** Curated survey counts are inconsistent across the main text, the survey file, and the appendix: 89, 92, and 93.
- **M3.** The strong-regularization bootstrap coverage is quoted as 0.338; the supplied summary gives 0.339.
- **M4.** The large-bias shortcut in the critical-value function returns bias plus the one-sided quantile times the standard deviation, which is anti-conservative by a negligible amount. Say so or remove the shortcut.
- **M5.** The neural checkpoints selected are the last available epochs on two stacks, so the neural baseline may be undertrained. Irrelevant to uncertainty, but should be stated if the study is kept.
- **M6.** Report the shift budget in ångströms in the main text so readers can see it is negligible.
- **M7.** The "L2" adversary in figure captions needs a one-sentence definition in the caption itself.
- **M8.** Cite Low's 1997 fixed-length interval result and Cai and Low on adaptive confidence intervals when discussing the 1.2 factor and the impossibility of adaptation.
- **M9.** The related-work appendix is a long survey largely disconnected from the results and should be cut to what the experiments use.
- **M10.** The protocol's "retain every case" conflicts with the evaluator's abort-on-undercoverage assertion; record failures instead of raising.
- **M11.** State plainly in the abstract that all uncertainty experiments use simulated noise and known simulation limits.
- **M12.** The main text is far longer than an ICML main body; much of Sections 5 and 6 belongs in the appendix.

**Revision plan.** Essential items are numbered; optional items are lettered. Do not add further synthetic generators, angles, or stacks; the existing 144 settings already saturate what simulation can show.

1. **Reframe and cut.** Address R5 and R6 and M1, M9, M11, M12. Lead with the dictionary-exclusion demonstration and the continuous audit; drop record counts and the 1.2 factor from the abstract; move the prediction study out.
2. **Pose-aware continuous optimization.** Address R3 on the existing three geometries at 0.5°, 1°, and 2°, with the gauge decision made explicit. Report widths, feasible-over-upper ratios, and sign power.
3. **One end-to-end experimental attempt.** Address R1 on EMPIAR-10028, with noise covariance from held-out exposures, an externally sourced pose radius, a declared density prior, and a realistic translation budget. Report the result whatever it is.
4. **Scale test.** Address R4 with a matrix-free pose audit at ten thousand particles and a higher band.
5. **Realistic nuisance model.** Address R8 within the experiment from item 3.

- **a.** Ablate a positivity or support-constrained class to quantify how much width the signed radius-two class costs.
- **b.** Add a descriptive half-map or local-resolution figure for the same feature, labeled as a different estimand, to give readers a familiar reference point.
- **c.** Drop the Monte Carlo cross-checks and binomial intervals; they add nothing to analytic coverage.

**Verdict.** Reject. Confidence 4 of 5. Is the current work a strong contender for ICML acceptance? No. The reasons are, in order of weight: the method cannot be applied to experimental data with its guarantee intact and the paper's own diagnostics show why; the demonstrated intervals are informative only for very coarse features under pose budgets tighter than real alignment achieves; the one new bound is loose by a factor of two to six at realistic budgets and is applied to weights not optimized for it; and the statistical contribution is a specialization the authors themselves describe as classical. What is good, and should be preserved in any resubmission, is the analytic continuous audit itself, the careful proofs, the dictionary-exclusion demonstration, and the reproducibility discipline. If items 2 and 3 of the plan produce nonvacuous experimental intervals, this becomes a credible paper; if they do not, it becomes an honest negative result that would need a different venue and framing.
