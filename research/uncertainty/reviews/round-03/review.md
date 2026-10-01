I'll start by inspecting the review directory and the evidence copies before forming any judgment.

# Independent review, round 3 — *Continuous Audits of Representation and Alignment Uncertainty in Cryo-EM* (focused manuscript, 26 pp.)

**Verdict: reject. Confidence 4 of 5. Strong contender for ICML acceptance: no.**

The mathematics I checked is correct, the code I read implements it, and every number I traced matches its record. The manuscript is candid that its pose-aware intervals are vacuous. The rejection rests on what remains:

- The fixed-pose audit is a ridge estimator with a bias-aware critical value. The matched Gaussian-prior interval shares its guarantee and is slightly narrower on two of three stacks.
- The pose-aware bound is uninformative for structural reasons that the manuscript misdiagnoses.
- The 600-dataset study cannot discriminate between procedures, and its outcome for the proposed procedures was fixed by its design.
- The one finding that study does contain, a clear estimated-pose bias, is not reported.

A model review does not predict a conference decision.

## 0. What I examined

**Pages.** All 26 rendered pages were visible: main text 1–6, references 7–10, appendices 11–26. I re-opened pages 6 and 14 from `pages/`.

**Method.** Reading only. I executed nothing and browsed nothing.

| Category | Read |
|---|---|
| Code | `run_uq_end_to_end_local.py`, `calibrate_uq_local_alignment.py`, `summarize_uq_end_to_end_local.py`, `uq_end_to_end(_summary).py`, `uq_mixed_pose.py`, `uq_local_alignment.py`, `uq_continuous.py`, `uq_continuous_quadrature.py`, `uq_continuous_gaussian.py`, `uq_cached_quadrature.py`, `uq_intervals.py`, `uq_fourier_pose_baseline.py`, `benchmark_uq_continuous_gaussian_v2.py` and its summarizer, `apply_uq_fresh_noise.py`, `uq_noise_calibration.py`, `uq_group_noise.py`, `uq_centered_noise.py`, `report_uq_breakdown_radii.py`, `uq_joint_bias.py`, `uq_phase_control.py`, `run_uq_phase_control.py`, `uq_cell_moments.py`, `uq_continuous_pose.py`, `uq_data.py`, `audit_uq_grid_refinement.py`, `benchmark_uq_fourier_pose.py`; `uq_cubic_pose.py` skimmed |
| Tests (read, not run) | mixed pose, local alignment, continuous Gaussian, end-to-end, phase control |
| Refitting records | All three summaries (error distributions for all 72 cells per stack, alignment and convergence blocks); replicates 10049/000, 10028/000, 10028/137, 10076/117; pattern scans over all 600 replicate files |
| Other records | Calibration summaries and two replicate records; Gaussian v2 files (all 48 class-coverage values, two features' 14 checks each); the interrupted v1 record; all 48 fresh-noise rows; centered calibration for 10076; breakdown summary; registered replay for 10028; phase-control summary (part); RELION 10049 metrics; frozen prediction metrics; the 10,000-particle audit |

**Not verified.**
- Array files (weights, poses, pre-iteration objectives) were not available to me.
- I could not check git ordering of protocol freezes, or the reported test-suite status (279 passed, 1 skipped).
- I read no third-party papers. Many 2026 bibliography entries are unverifiable by me.
- Figures 2–6 were checked only visually and against their JSON summaries.
- The two older frozen studies were checked at summary level only; the current manuscript no longer reports them.

## 1. Contribution and closest prior art

**What is delivered.**

1. **Fixed-pose continuous audit.** The worst-case bias of an affine estimator over an L² ball on the unit cube, computed from a sinc-kernel Gram and truncated-Gaussian transforms, with a quadrature remainder and a feasible dual.
2. **Normal-envelope proposition.** For z ≥ √3, the interval z√(s²+b²) covers under any fixed bias up to b. A Gaussian prior with scale τ = B therefore shares the fixed-pose uniform guarantee. A companion counterexample covers noise-adaptive bias.
3. **Conditional mixed-pose bound.** Hoeffding on first-order pose terms plus an integrated second-order remainder, with a maximum-score calibration and a one-frequency phase model.
4. **Negative empirical results.** A 600-dataset local-refinement study, a 48-fit matched Gaussian comparison, experimental breakdown radii, and a registered dictionary replay.

**Closest supplied prior art.** Donoho (1994) and Armstrong–Kolesár (2018) supply the construction, critical value and duality. Backus–Gilbert (1968) supplies resolving kernels. Ullrich et al. (2020) and Rangan et al. (2024) supply Fourier-slice density uncertainty and pose–density coupling. Jensen (2001) and Shaikh et al. (2003) supply alignment attenuation and cross-validation; the authors read abstracts only. The authors credit these and claim no new principle.

**Prior art not in the packet (my own background, not verified here).**
- The τ = B result is the known link between Gaussian-process posterior variance and worst-case error over a norm ball (Kanagawa et al. 2018; Srinivas et al. 2010; Fiedler et al. 2021).
- Section 5 and Appendix D describe the incidental-parameter problem (Neyman–Scott 1948), as studied in multi-reference alignment (Perry et al. 2019) and as reference bias (Henderson 2013).
- An ICML audience will expect both connections.

**Net novelty.** Operator-specific integrals and a sharp elementary constant. There is no method that outperforms a baseline in any tested regime.

## 2. Technical checks

I re-derived each item by hand and compared it with the code. All are correct as stated.

| Item | Note |
|---|---|
| Theorem 1 | Premise now explicitly excludes noise-dependent nuisances |
| Proposition 1 and sharpness | Termwise series argument holds; I reproduced the coefficient φ(z)z(z²−3)/6 |
| Adaptive-bias counterexample | Infimum 2Φ(√(c²−1))−1; 0.908 at c = 1.96 |
| Duality, feasible point, solve-gap identity, preconditioner ordering | Weak duality and P ⪯ H verified |
| Sinc Gram and truncated-Gaussian transform | Factor ½, sum-frequency term and complex-argument CDF all correct |
| Quadrature remainder | Reproduced 3.2×10⁻²³ at radius 5, order 40 |
| Cell moments | Series coefficients in `normalized_cell_moments` checked term by term |
| Pose derivative Gram; moments 1/12 and 13/720; remainder (23) | 13/720 attained at b = (1,1,1)/√3 |
| Theorem 2 | Correct under its premise |
| Maximum-score calibration | 0.04 + 1/129 < 0.05 holds marginally |
| Phase model, generating function (24), bound (25) | Re-derived through the chi-square decomposition |
| Noise lemmas (trace, group, Helmert centering) | Correct; premises unverified |

Four findings change how the results should be read.

**A. The fixed-pose audit is a ridge estimator.**
- Stationarity of z‖w‖ + B‖ℓ − A*w‖ gives (G + λI)w = a with λ = zb/(B²s). The audit weights are therefore Gaussian-prior weights at τ_eff = B√(s/(zb)).
- From the records, τ_eff ≈ 0.41 (10028, radius 5), 0.81 (10028, radius 12) and 3.4 (10049, radius 12).
- "Audit versus Gaussian baseline" compares two points on one ridge path with two critical values.
- On 10049 and 10076 the τ = B point is narrower even with the cruder critical value (ratios 0.994 and 0.999). The sum-objective weights are thus measurably suboptimal for the reported width.
- Fix: minimise the folded-normal width over τ in one dimension and report that.

**B. The pose bound is vacuous for structural reasons, not because of the second-order term.** The manuscript attributes 99.3–99.5% of the width to the nonlinear remainder. That arithmetic is right; the implied diagnosis is not.
- **First-order term alone.** Across all 2,400 fits, the smallest Hoeffding term is 48.8 against a no-data half-width of 16.2 (10076, replicate 106). With the remainder set to zero, every procedure would still fall back.
- **Any per-frequency bound.** In 10049 replicate 0 the remainder is 16,133. That implies Σ|c| ≥ 150, so the crudest valid bound, (B+‖ρ₀‖)·Σ|c|·min(2, rL₂) with no Taylor expansion, is at least 628. That is 39 times no-data.
- **Scaling.** The remainder scales with radius squared. It reaches half of no-data only if the worst per-particle error is below about 0.7–0.9° and 0.35–0.45 Å. That is tighter than the study's own 2° initialization.
- **Realized effects.**
  - In the radius-12 Gaussian records, a 2° perturbation of every particle moves the expected centre by 0.005 on a feature of 1.63 (10049, region 1) and by 0.027 on 3.76 (10028, region 1).
  - The certified one-degree bias exceeds the no-data half-width of 35.5 on 10028.
  - The gap is about three orders of magnitude.
- **Attenuation.** With centred errors, E[e^{iφ(ξ)}] ≠ e^{iφ(0)}. Theorem 2 charges this to the remainder in the worst case. An audit built on the mean operator E_ξ A(ξ) would have no remainder, at the price of declaring the error law.
- **Cheap fix that does not rescue it.** Replace the remainder per frequency by min{r²(L₄²+H₂)/2, rL₂ + min(2, rL₂)}. This is about 2.7 times tighter at the calibrated radius and still hundreds of times no-data.

**C. Theorem 2 was never exercised where its premise holds.**
- The remark that E[u | θ̂] = u is correct for fixed true poses. It is the reason the standard model treats poses as random.
- In that model, with a fixed template, posterior-mean pose errors are centred and independent across particles given the alignment images.
- The study's generator has this structure: truth is Gaussian about the supplied initialization.
- The estimator ignores that structure and runs unregularized maximum-likelihood refinement.

**D. Calibration is marginal, and one batch shows it.**
- On 10028 with the oracle template, 7 of 200 test scores exceed the calibration maximum (inclusion 0.965).
- Under exchangeability the chance of seven or more is ∏ᵢ₌₀⁶(200−i)/(328−i) ≈ 3.0%. That is unremarkable across six cells.
- Conditional on that batch, the nominal error is about 0.04 + 0.035, which exceeds 0.05.
- A batch-conditional guarantee at 1% with 95% confidence needs about 298 calibration datasets, not 128.

## 3. Experimental audit

### 3.1 Three kinds of evidence

| Level | What exists | What it shows |
|---|---|---|
| Conditional simulation coverage | 672 analytic Gaussian checks at fixed weights and prescribed poses; 12 dictionary replays | Implementation verification, and one illustration that a dictionary-restricted interval ignores bias |
| End-to-end simulation | 600 datasets, 216 cells, all 200/200 | Nothing about validity; see 3.2 |
| Experimental prediction | Held-out image error; CryoLike scores | Unrelated to density coverage, as the authors say |
| Experimental density coverage | None | 48 intervals under assumed inputs; five disagreements with the registered 10076 reference under centered calibration |

### 3.2 The 600-dataset study

Round 2 specified radius 5 and 128 particles, and the authors ran that faithfully. The records show the specification was inadequate as a test bed.

**The alignment step fits noise.**
- Each particle supplies 80 real coordinates. Mean squared whitened signal is 1.02, 0.13 and 0.11 per coordinate on the three stacks.
- In 10049 replicate 0 (oracle template) the objective falls at every iteration, from 10,137 after iteration 1 to 9,638 after iteration 10. Pose error grows from 2° to 15.5° RMS.
- The final value is about 600 below the noise-only expectation of 10,240, roughly five per particle for five fitted parameters. The pre-iteration value is in arrays I could not open.
- Across all 12,000 iteration records, only two show fewer than 110 of 128 particles still accepting a step. The fits are truncated, not converged.
- No particle reaches the declared 20° or 10 Å search limits. Travel is set by ten iterations of 2° and 1 Å steps.
- In one calibration record the rotation-error median is 16.2° and the 90th percentile is 20.1°.
- The "calibrated" radius of about 30 therefore measures the optimizer's caps.

**The outcome for the pose procedures was predetermined.**
- If poses were simply left at the initialization, the calibrated score would be about 5.5.
- The remainder would then be about 500 on 10049, which is 31 times no-data.
- No data needed to be generated to learn that these procedures always fall back.

**The coverage endpoint has no power.** The table covers the fixed-pose audit estimator. The four mean errors are for oracle template same image, oracle independent, pilot template same, pilot independent.

| Stack, target | Truth | Pilot value | Half-width | Noise SD | Mean error (O-same, O-indep, P-same, P-indep) | Largest error ÷ half-width |
|---|---|---|---|---|---|---|
| 10028 centre | 2.73 | 0.88 | 1.40 | 0.10 | −0.03, −0.06, −0.12, −0.13 | 0.31 |
| 10028 contrast | −0.61 | −0.16 | 2.06 | 0.14 | −0.02, +0.02, +0.03, +0.03 | 0.22 |
| 10049 centre | 3.04 | 2.81 | 1.82 | 0.30 | +0.21, −0.09, +0.19, −0.17 | 0.54 |
| 10049 contrast | 0.50 | 0.86 | 2.51 | 0.29 | −0.04, +0.20, +0.30, +0.41 | 0.53 |
| 10076 centre | 3.58 | 4.41 | 1.72 | 0.25 | +0.09, −0.10, +0.24, −0.08 | 0.59 |
| 10076 contrast | 0.21 | 0.60 | 2.45 | 0.31 | +0.05, −0.01, −0.05, +0.03 | 0.35 |

- Half-widths are 6 to 14 noise SDs because the B = 2 density budget dominates.
- The worst error in 200 trials never exceeds 59% of the audit half-width. For the narrowest comparator (Gaussian, τ = B/2) the worst is 87%.
- The 200/200 result reflects the density budget absorbing everything on a benign generator.
- The authors' own earlier frozen study put fixed-pose coverage as low as 0.67 on in-class boundary densities under 1–2° coherent errors. The present study has no such arm.

**What the data do show, and the manuscript does not report.**
- **Estimated-pose bias is detectable.** Mean errors reach 1.2–1.4 noise SDs (10028 centre with pilot template; 10049 contrast with pilot template). The Monte Carlo standard error is about 0.07 SD.
- **Noise-only coverage would fail.** By my normal approximation, a noise-only 95% interval would cover about 71–75% in those cells.
- **Same-image selection.** Same-image centres exceed independent-image centres by 0.30 on 10049 and 0.19 on 10076, about ten standard errors. This is the paper's phase-model mechanism appearing in its cryo-EM simulation.
- **Reference bias.** On the 10049 contrast, pilot-aligned estimates average 0.80–0.91. That sits at the pilot's own value of 0.86, not at the truth of 0.50.

**Missing arms.** There is no true-pose control, no noise-only comparator and no in-class adversarial generator. Without them the study cannot separate representation bias, attenuation and noise selection.

### 3.3 Frozen versus development

- The runners refuse to start unless protocol and sources match HEAD, and they save failures.
- The protocol states that calibration pose errors were inspected before the interval rules were frozen. The degenerate alignment was therefore known in advance.
- The remainder decomposition, breakdown radii and Gaussian rank amendment are post-outcome and labelled so.
- The interrupted rank-1,024 run is preserved: four fits, none converged.
- Table 9 shows development-group prediction errors although a frozen fresh-exposure evaluation exists. Its 10028 values are 0.8353, 0.8414 and 0.8294, against 0.8355, 0.8417 and 0.8293 in the table.

### 3.4 Baseline fairness

- The matched Gaussian comparison is fair in operator and scale.
- The record holds worst-case class coverage for τ = B/2: 0.980 on 10028, 0.9949 on 10049 and 0.9954 on 10076, against nominal 0.99583. The manuscript reports only generator-specific minima.
- At τ = B all 24 values are at least 0.9958364, as Proposition 1 requires.
- The first-order Gaussian pose term adds 0.03% to variance on 10049 at radius 12. It is not validated either.

### 3.5 Calibration

- The experimental noise lemmas are sound under common Gaussian covariance and independent exposures. Nothing in the packet supports those premises.
- The noise bound on 10028 is 0.28, against the 0.15 the weights were designed for.
- The single raw movie shows neighbour correlations near 0.66 and establishes no noise law.
- The six, four and zero exclusion counts, and ten and six under centering, match the records.

### 3.6 Usefulness and scope

- The only sign any procedure resolves in the 600 datasets is that of the central average of a centred particle, which is known beforehand.
- Contrast sign exclusion is at most 0.005.
- On 10049 the pilot alone is within 0.22 of the truth. That is less than the data-based estimator's noise SD of 0.30.
- Reported widths of 0.09–0.13 of no-data are ±50–60% of the central value and ±340–1,200% of the contrast.
- At radius 12 the bias bound on 10049 is 0.027 against a noise SD of 0.225. The representation audit adds about 0.6% to a noise-only interval there.
- No interval excludes zero under any rotation budget.
- Three acquisition designs support no generalization claim.

### 3.7 Statements not supported as written

1. **"Large nonlinear remainders dominate."** This implies the remainder is the obstacle. Finding B shows it is not.
2. **"Conditional centered-error bounds can be informative under their model."** No regime in the packet shows this. By the scaling above it needs worst-case errors below 1°.
3. **"Re-estimates local poses."** The records show truncated noise-driven drift bounded by step caps.
4. **"Exposes conservatism rather than a coverage advantage."** It exposes conservatism of the density budget on benign generators. It cannot speak to estimated-pose coverage within the class.
5. **The alignment thesis.** Within the manuscript it rests only on a one-frequency toy model. The supporting cryo-EM evidence sits unreported in the authors' own summaries.

## 4. Concerns

### Status of round-2 concerns

| ID | Concern | Status |
|---|---|---|
| R1 | Experimental inputs assumed | Open; one movie analysed, nothing validated |
| R2 | No useful inference | Open |
| R3 | Loose pose certificate | Superseded by R14 |
| R7 | Novelty | Open; conceded in text; sharpened by finding A |
| R8 | CTF and scale nuisances | Open; not addressed |
| R9 | Wrong pose model | Lemma added; premise never met; criterion failed |
| R10 | Matched baseline | Done; outcome is near-equivalence |
| R11 | End-to-end coverage | Run as specified; criterion failed (642–1,242 against 0.5); study uninformative |
| R12 | Reporting | Largely fixed; refitting tables still omit widths relative to the estimate |

### Major concerns

**R11. No informative end-to-end coverage evidence.**
- *Evidence.* Section 3.2.
- *Change.* One redesigned simulation with:
  - a regime where a regularized or posterior pose estimate improves on initialization, checked by a pre-study gate;
  - a true-pose arm, a noise-only arm and an in-class boundary generator;
  - at least 200 datasets on three geometries.
- *Acceptance.* The gate passes; at least one comparator visibly undercovers; the proposed interval covers at nominal within Monte Carlo error.
- *Status.* Fixable. Essential.

**R13. The study's real finding is unreported.**
- *Evidence.* Biases of 0.3–1.4 noise SDs; same-versus-independent centre shifts; reference bias.
- *Change.* From stored arrays, report standardized bias, RMSE against the pilot, and coverage of noise-only and reduced-B intervals. Add a true-pose arm.
- *Acceptance.* These appear in the main text, with one figure.
- *Status.* Fixable in days without new simulation.

**R14. The pose-aware bound is structurally vacuous and misdiagnosed.**
- *Evidence.* Finding B.
- *Change.* Take one of two routes:
  - prove a lower bound showing the class forces near-no-data widths at a given pose radius; or
  - replace per-frequency triangle bounds and the L² ball. Compute the exact class-worst-case bias at realized poses with the existing sinc Gram as a diagnostic, then build on a mean operator or a smaller physical class.
- *Acceptance.* A theorem, or an interval below 0.5 of no-data and within twice the fixed-pose width at errors a working estimator attains.
- *Status.* Fatal to the pose-aware contribution as formulated. Fixable only by a different construction.

**R9. Theorem 2's premise is never satisfied or tested.**
- *Evidence.* Finding C.
- *Change.* Adopt the latent-pose model explicitly. Use posterior-mean poses or marginalization with a declared prior. Test Theorem 2 where it applies.
- *Acceptance.* Nominal coverage in the R11 study with the premise holding by construction.
- *Status.* Potentially fixable.

**R7. Novelty.**
- *Evidence.* Finding A; Table 1.
- *Change.* State the ridge-path equivalence. Run the one-dimensional τ search. Cite the kernel and incidental-parameter literatures.
- *Acceptance.* One pre-declared regime where the audit differs materially from the matched baseline and is shown to be the more trustworthy.
- *Status.* Not fixable by wording.

**R2. Usefulness.**
- *Evidence.* Section 3.6.
- *Change.* Demonstrate a non-trivial feature: a contrast, or a feature of SD ≤ 20 Å at radius ≥ 12.
- *Acceptance.* Sign power at least 0.8 on two stacks with nominal coverage in R11.
- *Status.* Fatal to a usefulness claim. Depends on R14.

**R15. Density class and metric.**
- *Evidence.*
  - The B = 2 ball is the trivial distance between unit-norm functions.
  - No-data widths are 9–70 times the pilot's actual error.
  - The class's adversarial densities drive every pose term.
- *Change.* Use a physically argued class (support mask, spectral decay, positivity) with a radius from independent data. Report widths relative to the estimate.
- *Acceptance.* Pose terms within twice fixed-pose at realistic radii under the new class.
- *Status.* Uncertain. Fixable in principle.

**R1. Experimental inputs.**
- *Status.* Unchanged from round 2. Essential only if experimental claims return.

### Minor concerns

- **M1 (carried).** Not anonymized: repository URL and project name appear on page 1.
- **M8 (carried).** Low (1997) is uncited. Cai–Low (2004) is in the bibliography file but no longer cited.
- **M23.** Add worst-case class coverage for τ = B/2 to Table 1.
- **M24.** Table 9 should show the frozen evaluation, or be labelled as development only.
- **M25.** The mixed interval does not reduce to Theorem 1 as radii go to zero: the multiplier is 2.72 rather than a normal quantile. Split α between a Gaussian quantile and the Hoeffding term.
- **M26.** State the batch-conditional calibration failure on 10028.
- **M27.** The reading note for "Beckers, Palmer and Sachse" points to a 2020 article. The cited entry is Beckers, Jakobi and Sachse (2019).
- **M28.** Layout:
  - no figure in the main text;
  - page 6 is about 60% blank;
  - tick labels overlap in the lower row of Figure 1;
  - Tables 4–6 spend three pages on entries of 1.000;
  - Appendices E.6, E.7 and F are tangential.
- **M29.** Local-alignment tests cover only the noise-free case.
- **M30.** Take the minimum of the remainder with the elementary bound in finding B.

## 5. Revision plan

**Essential, in order.**

1. **Re-analysis of existing arrays (R13).** Add the exact class-worst-case bias at realized poses. This decides whether vacuity is intrinsic or an artefact of the bound.
2. **Reformulate or prove (R14, R9).** Either a lower bound that turns the negative result into a theorem, or a latent-pose construction that is non-vacuous.
3. **One redesigned end-to-end simulation (R11).** Include the gate and the three missing arms.
4. **Rewrite (R7).** Cover the ridge-path equivalence, the τ search and the missing literature. Put one main-text figure of the bias finding in place of the 1.000 tables.

**Conditional.** One homogeneous experimental stack with evidence-based inputs (R1), only if experimental claims are kept.

**Optional.** Higher-order remainders; CTF and scale nuisances; optimized weights at 10,000 particles; more geometries.

**Not needed.** More radius-5 replicates of the current design; further calibration variants; further reconstruction baselines; survey growth.

## 6. Verdict

**Reject. Confidence 4 of 5.** I read the theory, the code paths behind every headline number and a broad sample of raw records. I ran nothing and could not open arrays.

**Is the current work a strong contender for ICML acceptance? No.** Reasons, in order of weight:

1. No procedure proposed here is informative at any pose-error level the evidence supports, and finding B shows why refinement of the same bound cannot change that.
2. The fixed-pose audit is matched or beaten by a Gaussian baseline with the same guarantee.
3. The flagship study has no power and a predetermined primary outcome.
4. The experimental inputs remain unvalidated.
5. The statistical content is classical.

**Worth preserving.** The sharp envelope proposition with its adaptive-bias counterexample, the honest matched baseline, the failure-retaining runners, and the provenance discipline.

**What would change the assessment.** A lower bound, or a latent-pose construction that is validated and non-vacuous in a regime where alignment works. Reporting the bias finding already in the data would, by itself, make the manuscript a more useful negative result.
