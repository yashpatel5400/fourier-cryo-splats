Re-checking the manuscript text and tables against the raw records one last time so every number in the review is cited accurately; the written review follows.

Spot-checking the centered-calibration and baseline records, then writing up the review.

Verifying one last RELION figure for 10028, then I'll deliver the review.

Nearly done: confirming two amplitude values, then writing the review itself.

Reading the cubic and pose-optimization modules I hadn't yet covered so the implementation assessment is complete.

# Independent review, round 2 — *Continuous Fourier-Slice Audits for Cryo-EM Density Uncertainty* (v0.6.2-dev, 55 pp.)

**Verdict: reject. Confidence: 4 of 5. Strong contender for ICML acceptance: no.**

The mathematics I checked is correct, the code I read implements it, and every number I traced to a record matches. The rejection rests on what the work delivers:

- At fixed poses, the method is a closed-form continuous version of classical ridge / optimal-recovery inference. Its intervals sit within tens of percent of a properly matched Gaussian-prior interval.
- With rotation budgets of 1° or 2°, none of the 24 experimental intervals excludes zero, and 16 revert to the no-data interval.
- There is no end-to-end coverage evidence, not even in simulation.
- The experimental inputs remain assumptions, and the project's own diagnostics contradict some of them.

The disclosures are careful, but a disclosed limitation is still a limitation.

## 0. What I examined and what I could not verify

Everything below comes from reading. I executed no code, ran no tests and browsed nothing.

| Category | Read |
|---|---|
| Manuscript | All TeX sources (`main.tex`, every `\input` file and table); rendered appendix pages 20–55 |
| Core code | `uq_intervals`, `uq_continuous`, `uq_continuous_quadrature`, `uq_continuous_pose`, `uq_continuous_moments`, `uq_pose_operator`, `uq_pose_optimization`, `uq_pose_exchange`, `uq_shift_only`, `uq_cell_pose_pair`, `uq_random_spectral`, `uq_joint_bias`, `uq_ball_remainder`, `uq_higher_remainder`, `uq_cubic_pose`, `uq_cubic_design`, `uq_cubic_optimization`, `uq_cubic_subspace`, `uq_cubic_enrichment`, `uq_joint_cubic_design`, `uq_trust_region`, `uq_two_pose_modulus`, `uq_noise_calibration`, `uq_group_noise`, `uq_centered_noise`, `uq_fourier_variational`, `uq_fourier_pose_baseline`, `uq_data`, `physics` |
| Scripts | `apply_uq_fresh_noise`, `apply_uq_pilot_targets`, `run_uq_pilot_targets`, `confirm_uq_continuous`, `audit_uq_grid_refinement`, `benchmark_uq_fourier_pose`, `probe_uq_experimental_noise`, `diagnose_uq_background` |
| Tests (read, not run) | continuous, continuous_pose, cubic_pose, pose_optimization, fourier_pose_baseline, group_noise, noise_calibration, noise_cohort |
| Raw records | Fresh calibration (summary and all three stack files); centered calibration (10028, 10049 in full; 10076 by search); two locked fixed-pose fits; one enclosing-domain audit; one cubic-weight probe; the 10,000-particle audit; 10028 reference registration; 10028 RELION metrics; background diagnostics; pilot and reference amplitudes |

Not verified:

- **Execution.** The "246 passed, one skipped" test status is the authors' report.
- **Frozen studies.** For continuous-v1 and continuous-moments-v2 the packet holds only summary files, not per-setting records.
- **Figures.** No figure files are in the evidence directory, and rendered pages 1–19 were not visible to me, so I have not seen the main-text figures.
- **Freeze ordering.** I have no git access to confirm commits 76dbd53, 8d66785, 60efd9b or becc3a6 preceded their outcomes.
- **Third-party papers.** I read none. Statements about them rest on my background knowledge and the project's annotations. I cannot confirm many of the 2026 bibliography entries.
- **Unread code.** The mixture-likelihood modules, CryoLike adapter, RELION and cryoDRGN pipelines, earlier finite-voxel modules, `uq_cubic_coordinate_solver` and `uq_cubic_design_api`.
- **Omitted data.** Optimization traces and weight vectors.
- **Floating point.** NUFFT and roundoff accuracy.

## 1. Contribution and closest prior art

**What is delivered.**

1. **Fixed-pose continuous audit.** For an affine estimator centred on a pilot, the worst-case bias over an L² ball on the unit cube is B‖ℓ − A*w‖. It is computed exactly from a sinc-kernel Gram and truncated-Gaussian transforms. A matrix-free Gauss–Legendre variant carries an explicit remainder and a primal–dual bracket.
2. **Pose audit.** Each particle has a five-dimensional pose ball. The audit uses a shared-field spectral bound for first- and second-order fields, integrated-moment and enclosing-domain remainders, a full cubic lift, cross-term and known-pilot refinements, a randomized spectral certificate, and pose-aware weight optimization with a dual bound that does not depend on optimizer success.
3. **Noise-scale lemmas.** A noncentral trace bound and its directional, group-inflated and centered versions.
4. **Diagnostics.** Two-pose modulus lower bounds, feasible adversaries, and CTF, support and sign ablations.
5. **Experiments.** Four frozen protocols and a large development record on three EMPIAR geometries, plus matched Gaussian-prior baselines and RELION, CryoLike and cryoDRGN context.

**Closest supplied prior art.**

- Donoho (1994) and Armstrong–Kolesár (2018) supply the statistical construction, the duality and the efficiency factor.
- Cai–Low (2004) supplies the two-point bounds.
- Ullrich et al. (2020) supply probabilistic Fourier-slice reconstruction and the variance-versus-bias distinction.
- Rangan et al. (2024) supply pose–density soft modes.
- Batlle et al. (2025) cover simultaneous functionals.
- Wasserman et al. (2020) and Lindsay (1983) underlie the mixture appendix.

The authors credit all of these. Two connections are not cited; both come from my own background, not the packet. The fixed-pose construction is the Backus–Gilbert resolving-kernel picture. Equivalently, it is Gaussian-process regression with a white L² prior on the cube (finding B below).

**Net novelty.** Operator-specific integrals, Taylor constants and software. The statistical principle is not new, as the manuscript says.

## 2. Technical checks

I re-derived each item by hand and compared it with the code.

| Item | Status | Note |
|---|---|---|
| Theorem 1 | Correct as stated | Needs nuisances that do not depend on inference noise; finding A |
| Dual and feasible point v = t·h | Correct | Weak duality verified; code matches |
| Sinc Gram, factor ½, sign of sum-frequency block | Correct | |
| Truncated-Gaussian transform | Correct | Completing the square; complex normal CDF with consistent sign |
| Gauss–Legendre remainder | Correct | Reproduced 3.2e-23 (radius 5, order 40); about 5e-33 at radius 12, order 80, before coefficient scaling |
| Polynomial-field quadrature pad; √2 Gram-action factor | Correct | |
| Pose derivative fields | Correct | Phase Hessian a²[(k_a x_b + k_b x_a)/2 − δ_ab k·x]; symmetric vectorization preserves the unit ball |
| Shared-field spectral inequality | Correct | Elementary |
| Cubic remainder L³ + 3LH + T; moment constants | Correct | Uniform moments 1/12, 1/80, 1/448 |
| Sharp cube moments 13/720 and 205/36288 | Correct | Recomputed both |
| Cross-term, known-pilot and product-set bounds | Correct | |
| Enclosing-domain remainder; ball and cube kernels | Correct | (−1)^m sign verified |
| Cubic lift (55 columns), E₃ and Φ₃ | Correct | Code coefficients match term by term |
| Bell-polynomial higher-order remainder | Correct as a remainder | Not an interval, as stated |
| Randomized power bound | Correct | Predicts about 1.10 inflation at 4 probes and 40 steps; records show 1.09–1.16 |
| Pose-aware objective and dual certificate | Correct | Support argument survives log-sum-exp smoothing |
| Residual-controlled joint bound | Correct | |
| Two-point bounds, Corollary 1, two-pose modulus | Correct | Deterministic-length intervals only |
| Efficiency factor 1.198 | Correct at α = 0.05 | At the experimental α ≈ 0.00375 it is about 1.09 |
| Noise lemmas | Correct | The premise is the problem; finding D |
| Support isometry, tail bound, sign-class support function | Correct | |
| Gaussian posterior identities | Correct | Finding B |
| E-value and Lindsay bound | Correct | No usable test results |

Four findings change how the results should be read.

**A. The folded-normal step fails when a bounded nuisance depends on the inference noise.**
- Theorem 1 treats pose errors as fixed parameters. Consensus poses are estimated from the inference images.
- Counterexample: fix the design and weights, and let the pose error take whichever boundary value gives bias +b when wᵀε > 0 and −b otherwise. The error magnitude is then |wᵀε| + b.
- Coverage of ±q becomes 1 − 2α + 2Φ(−(q + b)/s), which tends to 1 − 2α when bias dominates. Miscoverage doubles even though the radius assumption holds.
- Fix for this submodel: use q = b + z₍1−α/2₎·s, which is the sum objective already being minimized.
- This does not repair the realistic case, where the design and weights themselves depend on the noise through the estimated poses. Only an end-to-end experiment can address that (R11).

**B. A matched Gaussian-prior interval carries the same uniform guarantee.**
- For α ≤ 0.05, the bias-aware critical value cv(t) at bias-to-noise ratio t = b/s satisfies cv(t) ≤ z₍1−α/2₎√(1 + t²).
- My check: a fourth-order expansion gives margin φ(z)·z(z² − 3)t⁴/6, positive when z² > 3. Direct evaluation at t = 0.2, 0.3, 0.5, 1 and 2 for α = 0.05 confirms it.
- With the exact continuous Gram and prior standard deviation τ = B per unit-norm direction, the posterior variance is ‖w‖² + B²‖ℓ − A*w‖². This is the manuscript's own identity.
- The credible interval is therefore uniformly valid over the B-ball. It is wider than the bias-aware interval by at most z√(1 + t²)/cv(t).
- That factor is about 1.0 where noise dominates (10049 and 10076 at fixed poses, where bias is under 2% of the half-width). It is about 1.4 at t ≈ 2 (10028).
- So "without the same uniform-class guarantee" describes the baseline as implemented, not Gaussian-prior inference.

**C. In the worst-case pose class, pose bias does not shrink with particle count.**
- Unbiasedness requires Σ A_i*w_i ≈ ℓ, so each ‖w_i‖ scales like 1/n and the noise SD like 1/√n.
- The worst-case pose terms are sums of per-particle norms of order 1/n, each of which an adversary aligns. Their total stays of order one.
- This is an analytical observation. The packet has no clean test of it: the 10,000-particle audit uses unoptimized weights, and its 0.655 width is dominated by a density residual of 6.54 against pose terms of 1.83 and 2.22.
- What the packet does show:
  - Independent Gaussian pose errors of 1° and 0.5 Å widen the matched baseline by 0.007–0.29%.
  - Two-pose lower bounds barely move with the pose budget: 0.087–0.147 at nominal poses, 0.087–0.151 at 2°, and 0.150 after the alternating probe.
  - Certified upper bounds at 2° are 0.43–0.63.
  - Feasible adversaries reach 0.27–0.36 of the bias bound.
- A factor of roughly three to seven separates what is known to be necessary from what is certified. The paper cannot say where the truth lies.

**D. The noise lemma bounds the calibration-average covariance, not the inference covariance.**
- Counterexample: if the inference exposures have covariance 2Σ₀ and the calibration exposures Σ₀, the true SD is 1.41 times the calibrated quantity.
- The lemma's built-in slack is only 1/√0.584 = 1.31 at n = 128.
- The project's own diagnostics show heterogeneity of this order. Per-image background SD at the 1st and 99th percentiles is 0.206 and 0.316 on 10076, a variance ratio of 2.3.
- The percentiles are 1/10/50/90/99, per `diagnose_uq_background.py:87`.

**Optimization and numerical claims.**
- The twelve locked fixed-pose fits are certified; the largest gap is 0.00437.
- The pose-aware and cubic designs are not. Recorded full-space gaps include 85%, 99.3% and 99.8% on pose-aware fits, and 0.927–0.977 on the cubic designs.
- The "restricted guide gaps" of 3.65e-4 and 8.85e-4 are not full-space gaps, as the text says.
- Integration error is bounded rigorously in real arithmetic. NUFFT tolerance, special functions and roundoff are checked empirically, not enclosed.
- A toy check of the quartic remainder found it 35.6 times the sampled error. This is disclosed.

## 3. Experimental audit

### 3.1 Three kinds of evidence

| Tier | What exists | What it shows |
|---|---|---|
| Conditional simulation coverage | Frozen v1 (96 settings), v2 (48 settings) and many development grids. Coverage is computed analytically from a deterministic bias and a known Gaussian SD, at prescribed poses inside the ball | Implementation verification. In-class coverage at or above nominal follows from Theorem 1; a failure would mean a bug. The informative outputs are width and sign power |
| Experimental prediction | Frozen held-out-exposure prediction error for Gaussian, voxel and cryoDRGN models | The neural model predicts slightly better on all three stacks (differences 0.0008–0.012 at NMSE 0.83–0.93). This has no bearing on density coverage |
| End-to-end density coverage | None | No run estimates poses from simulated images and then audits. No experimental truth labels exist |

The experimental unit is an acquisition design, so the effective sample is three stacks. It is two for a homogeneous-density claim, since 10076 is heterogeneous. Geometry subsets share exposures. Record counts such as 2,304, 960 and 384 are not replicates.

### 3.2 Frozen versus development

- Four protocols are frozen. The main text says "three" and then reports a fourth.
- In the fresh-exposure recalibration, only the SD bound was prospective. The code reuses the previously observed centre and bias bound (`apply_uq_fresh_noise.py:105-107`). The exclusion counts equal those of the earlier calibration.
- Everything else is development on inspected pools. This includes post-outcome choices: registration, broader priors, centered and projected calibration, and five cubic designs.
- The authors label these correctly and never take a minimum across alternatives.
- The frozen evaluator aborts on in-class undercoverage before saving (`confirm_uq_continuous.py:179-180`). The authors acknowledge this.

### 3.3 Baseline fairness

The Gaussian-prior baseline targets the same functional, so the comparison is legitimate.

- **Prior scale.** The quantity comparable to B is the prior SD per direction. The original priors have 0.0026–0.0106 per direction, 190–760 times tighter than B = 2. Their collapse to coverage 0.000 is uninformative.
- **"Much broader than radius two."** The RMS-norm comparison (19.0 and 189.6) depends on dimension. At coordinate SD 1 the prior is half of B per direction.
- **Discretization.** The baseline uses a 33³ trilinear grid with 16–25% forward discrepancy. The audit uses the exact operator. The same Gram could serve both.
- **Fixed-pose widths at coordinate SD 1:**

| Stack | Baseline | Audit |
|---|---|---|
| 10028 | 0.016–0.017 | 0.020 |
| 10049 | 0.062–0.074 | 0.053–0.054 |
| 10076 | 0.036–0.037 | 0.035 |

  "Narrower" holds on one stack at this prior. It holds on all three only at SD 0.1.
- **Under the tested 1° and 2° perturbations.** The baseline keeps fixed-reference coverage 0.986–0.9999 (nominal 0.9958) at those widths. At 1° the audit is at no-data on 10028, 0.46–0.47 on 10049 and 0.55–0.56 on 10076. At 2° it is at no-data everywhere.
- **Guarantee.** The audit never falls below nominal; the baseline does, modestly. The packet contains no realistic configuration where the baseline fails badly while the audit stays informative.
- **Worst cases.** Auditing the baseline's estimators over the signed L² ball gives widths of 0.257–0.815 of no-data, but those worst cases are signed functions that, as the authors say, need not resemble molecules.

### 3.4 Calibration

- **Cost of the lemma.** It inflates SD by 1.31, multiplies rows by up to √12 for exposure groups, and includes signal energy.
- **Mismatch with the design.** The resulting SD bounds are 1.9 times the simulated SD the weights were designed for on 10028 (0.283 against 0.151), and 4.55 times for the 10049 cubic estimator.
- **Premises.**
  - Background excess kurtosis is 0.18, 0.01 and 0.51 on the three stacks.
  - Lag-one pixel correlations of 0.20–0.29 are allowed by the directional lemma.
  - Heteroscedasticity across particles is not allowed (finding D).
  - The annulus may contain signal and ice, so these diagnostics cannot settle the question. Nothing in the packet supports the premise.
- **Registered references.** With centered calibration, five registered-reference values fall outside their intervals, all on 10076. This could be state mismatch or overconfidence; the data cannot tell.

### 3.5 Usefulness and scope

- **Feature scale.** Targets are Gaussian averages of SD 16–34 Å at radius 5, or 20 Å at radius 12, where the band ends at 19.7–40 Å.
- **Fixed-pose experimental intervals in feature units:**
  - 10028: 1.30 ± 1.07 twice (excluding zero), then 0.99 ± 1.06 and 1.01 ± 1.07.
  - 10076: 5.48 ± 3.04, 4.16 ± 2.84, 5.30 ± 3.04 and 4.93 ± 3.16.
  - 10049: none excludes zero.
- **The reporting metric.** The reported relative widths of 0.03–0.11 are ±55% to ±107% of the estimate. Excluding zero says a pilot-selected 20 Å blob has positive mean density.
- **With rotations.** No interval excludes zero.
  - At 1° the bias bound is 53–55 against a no-data half-width of 35.5 on 10028. It is 15.0–15.2 on 10076 estimates of 4.2–5.5.
  - At 2° the bias bounds are 33, 101–102 and 371–380 on 10049, 10076 and 10028.
- **Pose-error evidence in the packet.**
  - RELION's own estimates on the working data are 1.51°/2.0 Å, 8.87°/4.6 Å and 2.72°/2.9 Å.
  - Known-map information scales are 2.4°, 9.1° and 6.2°.
  - Deposited poses may be better, but the packet has no estimate of their error.
  - The 0.5° budget at which the simulated audit is informative is unsupported.
- **More particles.** With 1,024 particles and a 10 Å target, the fixed-pose relative width is 0.0138. Adding 1° and 0.5 Å gives 61.2 against a no-data 34.4. The best refinement reaches 0.751 with zero sign power.
- **CTF.** A ±100 Å defocus envelope sends 2 of 6 higher-band fits to no-data; ±500 Å sends all six.
- **Best pose-optimized design.** In simulation it reaches 0.175–0.180 of no-data, with full-space gap 0.93–0.98. On real images it gives 2.46 ± 3.01 to ± 4.30.
- **Density class.**
  - In experimental normalization, the recorded amplitude ratios put the reference norm at 0.26, 0.44 and 0.75.
  - The reference lies 1.01–1.04 from the unit-norm pilot (0.93–1.03 after registration).
  - The pilot is therefore no closer to the reference than the zero function is.
  - The class has no bandlimit, so B must also bound atomic-scale energy the data never see. No physical argument is given.

### 3.6 Claims not supported as written

1. **Abstract, "A frozen fresh-exposure noise recalibration gives intervals excluding zero for six of twelve".** The interval centres were observed before the freeze.
2. **Abstract, "narrower intervals, without the same uniform-class guarantee".** See finding B and section 3.3.
3. **"Much broader than a hard radius-two class."** Dimension-dependent.
4. **Main text, "All approximate-map values lie inside".** This uses the unregistered frame and raw calibration. The five registered disagreements appear only in the appendix.
5. **"Within a factor 1.20 of the minimum possible deterministic interval length."** This restates Donoho's bound plus solver convergence, at an α the experiments do not use.
6. **The dictionary-exclusion example (coverage 4.2e-5).** It was computed with the unregistered generator.

## 4. Concerns

### Status of round-1 concerns

| ID | Concern | Status |
|---|---|---|
| R1 | Experimental inputs are assumed | Open |
| R2 | No useful inference at relevant scale | Open; criterion not met |
| R3 | Pose bound loose; weights not pose-optimized | Partly addressed; criterion not met |
| R4 | Pose audit does not scale | Resolved as computation (10,000 particles, 4,965 s, 1.95 GB); says nothing about precision |
| R5 | Verification presented as calibration | Largely resolved; residue in R12 |
| R6 | Framing | Largely resolved; residue in M22 |
| R7 | Novelty | Open; sharpened by finding B |
| R8 | Nuisance scope | Partly addressed; results negative |

### Major concerns

**R1. Experimental inputs remain assumptions.**
- *Evidence.* Sections 3.4 and 3.5, and finding D.
- *Change needed.*
  - (a) Noise: either adopt a stationary noise-spectrum model with a stated validation check, or replace common covariance by a per-exposure scale model with exchangeability. Validate on a pre-declared split of the reserved cohort.
  - (b) Poses: estimate poses for inference images without their inference noise (frame-split or band-split), or drop experimental claims.
  - (c) Density: use a two-norm class, with a band-limited radius estimated from independent data and a physically argued total-energy bound. Report the breakdown radius at which each conclusion flips.
  - (d) Homogeneity: restrict experimental claims to 10028 and 10049.
- *Acceptance.* On one homogeneous stack, a non-fallback interval for a registered feature, with every input traced to evidence.
- *Status.* Fatal to any applied or experimental-coverage claim, and the main obstacle to significance. Potentially fixable with substantial work.

**R2. No useful inference at a relevant scale or realistic pose error.**
- *Evidence.* Section 3.5.
- *Change needed.* Demonstrate usefulness in the R11 simulation and on one experimental stack.
- *Acceptance.* Zero exclusion, or better a contrast between two regions, for a feature of SD ≤ 20 Å at radius ≥ 12, under a pose-error level supported by evidence. Sign power at least 0.8 on registered references in two stacks.
- *Status.* Fatal to a usefulness claim. Fixable only through R9 and R11.

**R3. The pose certificate is loose and the pose-aware designs are unconverged.**
- *Evidence.* Feasible-over-upper ratios are 0.27–0.36 at 2°, widths are 0.43–0.63, and gaps reach 99.8%.
- *Change needed.* If the worst-case class is kept, carry the Taylor polynomial to degree 4 or 5 inside the certificate (the remainder table shows 3.55 falling to 0.027). Bring full-space gaps on the 1° grid below 10%.
- *Acceptance.* The round-1 criterion: ratio at least 0.5 at 2° and widths below 0.5 on all stacks.
- *Status.* Fixable in principle. Optional if R9 is adopted.

**R7. Novelty.**
- *Evidence.* Finding B. The fixed-pose audit is continuous ridge or GP inference in closed form.
- *Change needed.* State the equivalences. Identify what the audit gives that a matched continuous Gaussian-prior calculation does not.
- *Acceptance.* At least one pre-declared regime where the audit's answer differs materially from the matched baseline and is shown, end to end, to be the more trustworthy.
- *Status.* Not fixable by wording. Contingent on R9–R11.

**R8. Nuisance scope.**
- *Evidence.* The CTF table. Per-particle scale is absent, and no per-particle radii exist.
- *Change needed.* Include CTF and scale errors in the R11 simulation, retaining cancellations.
- *Acceptance.* Non-vacuous radius-12 intervals at ±100 Å defocus error.
- *Status.* Fixable; lower priority.

**R9 (new). The deterministic worst-case per-particle pose class is the wrong nuisance model.**
- *Evidence.* Finding C.
- *Change needed.* Model alignment error as independent across particles and approximately centred, plus a small declared common-mode term.
  - First-order terms then enter as a bounded-difference variance, z·√(Σcᵢ²) in place of Σcᵢ, using field norms already computed.
  - Alternatively, constrain weights to be orthogonal to each particle's pose tangent and analyse estimated poses in a profile framework.
- *Acceptance.* On the twelve locked features at 1° RMS bounded error: half-widths within twice the fixed-pose values, sign power at least 0.8 on registered references in two stacks, and validity confirmed in R11.
- *Status.* Potentially fixable. This is the most valuable single change.

**R10 (new). The matched baseline is nearly equivalent at fixed poses and far more informative under tested pose perturbations.**
- *Evidence.* Finding B and section 3.3.
- *Change needed.* Compute the Gaussian-prior interval with the exact continuous Gram at τ = B and τ = B/2, with and without first-order pose marginalization. Run it on the locked features and in R11. Correct the text.
- *Acceptance.* The comparison is run with matched operator and scale, and the paper states the relationship plainly.
- *Status.* Fixable by one experiment and rewriting. The outcome may show the audit adds only a critical-value constant.

**R11 (new). No end-to-end coverage evidence.**
- *Evidence.* Section 3.1 and finding A.
- *Change needed.* One pre-registered simulation:
  - Radius 5, 128 particles, three geometries, two broad targets.
  - At least 200 independent replicate datasets per geometry.
  - Poses re-estimated from each simulated dataset by a fixed local refinement, and weights recomputed from the estimated poses. Fits take 6–14 s at this size.
  - The pose-error level set from a separate batch before coverage is examined.
  - Empirical coverage of the true feature reported for the fixed-pose audit, the pose-aware audit, the R9 interval and the R10 baseline.
- *Acceptance.* The proposed interval attains nominal coverage within Monte Carlo error at relative width below 0.5. Any shortfall or vacuity is reported.
- *Status.* Essential. Its absence is fatal to any end-to-end claim. Fixable in weeks.

**R12 (new). Reporting choices overstate informativeness.**
- *Evidence.* Section 3.6; the no-data-relative metric; unregistered references in the main figure.
- *Change needed.*
  - Report half-widths relative to the estimate and as breakdown radii.
  - Put registered references and the five disagreements in the main text.
  - State in the abstract that centres were previously observed.
  - Recompute the dictionary-exclusion example with the registered generator.
- *Acceptance.* The abstract and main text read correctly without the appendix.
- *Status.* Fixable by rewriting and one recomputation.

### Minor concerns

- **M1 (carried).** The manuscript is not anonymized. The repository URL identifies the author.
- **M2 (carried).** `SURVEY.md` calls 106, 113 and 119 the "current" count in different places.
- **M8 (carried).** Low (1997) on fixed-length intervals is still uncited. This is from my own knowledge.
- **M9, M12 (carried).** At 55 pages, the appendix is a development log.
- **M10 (carried).** The frozen evaluators still abort before saving a failure.
- **M13.** "Three completed protocols" should be four.
- **M14.** The 1.198 factor is quoted at α = 0.05. At the experimental α it is about 1.09.
- **M15.** The solver minimizes z‖w‖ + b but reports the folded-normal width. For the fixed-pose class a one-dimensional search along the ridge path would minimize the reported width directly.
- **M16.** Weights were designed for simulated noise and without group inflation. The noise-metric refit improved fixed-pose widths (0.111 to 0.098, 0.398 to 0.272, 0.158 to 0.133) and worsened 1° widths.
- **M17.** The 10028 pilot looks weak. The best registration correlation is 0.308, with a range of 0.175–0.308 over 48 starts. The locked-target panel looked banded rather than compact in the rendered page I saw.
- **M18.** Excess kurtosis of 0.51 on 10076 sits against the Gaussian premise.
- **M19.** Mark preprints, and tie each cited claim to a source that was actually read.
- **M20.** Qualify "certificate" wherever NUFFT and roundoff are unenclosed.
- **M21.** There is large float white space on appendix pages 47, 51–53 and 55.
- **M22.** The prediction, CryoLike, RELION and mixture sections are tangential to the uncertainty claim.

## 5. Revision plan

**Essential, in priority order.**

1. **End-to-end simulation (R11).** Include the matched baseline and the mixed-model interval. This single study separates conditional from end-to-end coverage.
2. **Mixed pose-error model (R9).** Build it from the existing field norms.
3. **Matched continuous Gaussian-prior baseline (R10).** Run it on the twelve locked features.
4. **One homogeneous experimental stack with evidence-based inputs (R1).** Report breakdown radii, whatever the outcome.
5. **Rewrite (R12, R7).** One thread, three evidence tiers stated plainly, development log moved out.

**Optional.**

- **a.** Higher-order certificate (R3).
- **b.** CTF and scale nuisances with cancellations (R8).
- **c.** Optimized weights at 10,000 particles.
- **d.** Positivity or support classes.

**Not needed.** Further generators, angles or stacks under analytic coverage. More calibration-contrast variants. More optimizer variants on the single cubic feature. Further mixture-likelihood bounding in this paper.

## 6. Verdict

**Reject. Confidence 4 of 5.** I read the mathematics and code closely and traced the numbers. I did not execute anything, did not see the main-text figures, and had no per-setting frozen records.

**Is the current work a strong contender for ICML acceptance? No.** Reasons, in order of weight:

1. No interval is informative at any pose-error level the packet's own evidence supports.
2. No end-to-end coverage has been tested.
3. A matched Gaussian-prior interval shares the fixed-pose guarantee and was far more informative under every tested pose perturbation.
4. The experimental inputs are unvalidated and partly contradicted by the project's own diagnostics.
5. The statistical contribution is a specialization of classical results.

**What should be preserved.** The analytic continuous audit, the real-arithmetic quadrature bounds, the dual certificates that do not depend on optimizer success, the registration audit that caught a real interpretation error, and the provenance discipline.

**What would change the assessment.** If items 1–4 of the plan yield validated, non-vacuous intervals under realistic alignment error, this becomes a credible submission. If they do not, it is a careful negative technical report that needs a different framing and venue.
