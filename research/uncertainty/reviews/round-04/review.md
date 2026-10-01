I'll start by inspecting the review directory and the exact evidence copies before forming any judgment.

# Independent review, round 4 — *Continuous Audits of Representation and Alignment Uncertainty in Cryo-EM* (33 pp., incl. Appendices G–J)

**Verdict: reject. Confidence 4 of 5. Strong contender for ICML acceptance: no.** A model review does not predict or guarantee a conference decision.

The mathematics I checked is correct, the code I read implements it, and every number I traced matches its record except one transcription slip. The rejection rests on what is delivered:

- **Main text.** Unchanged in substance since round 3. The fixed-pose audit is a ridge estimator with a classical critical value, the pose-aware intervals are vacuous, and the alignment-bias finding comes from an aligner that makes poses worse.
- **New branch (Appendices G–J).** It does not change the verdict. On two of three candidates the test cannot detect the 25% deletion it was designed for, at any sample size. It has no power once the viewing law may differ from uniform by a factor of two. It never touches an experimental image.
- **Coherence.** The main text never mentions the new branch, which studies a different estimand under the opposite stance toward the fitted map.
- **Direction.** A substantive change is needed. Section 5 gives a finite plan with a stopping rule.

## 0. What I examined

**Pages.** All 33 rendered pages were visible (main text 1–7, references 7–11, appendices 12–33). I re-opened pages 5, 32 and 33.

**Method.** Reading only. I executed nothing and browsed nothing.

| Category | Read this round |
|---|---|
| Core modules | `uq_moment_mc`, `uq_view_variance`, `uq_view_risk`, `uq_candidate_score`, `uq_fisher_score`, `uq_bispectrum`, `uq_preferred_views`, `uq_continuous`, `uq_refitting_diagnostics`, `uq_end_to_end` |
| Runners and verifiers | `probe_uq_candidate_scores`, `probe_uq_fisher_scores`, `probe_uq_replica_allocation`, `probe_uq_moment_mc`, `probe_uq_preferred_views`, `probe_uq_bispectrum`, `analyze_uq_candidate_scores`, `write_uq_candidate_score_results`, `analyze_uq_refitting_bias`, `verify_uq_fisher_scores`, `verify_uq_view_variance_sample` |
| Tests (read, not run) | candidate score, Fisher score, view variance, view risk, moment MC |
| Candidate study records | `variance-decomposition.csv` (all); scale-orthogonal paired-variance rows of `projections.csv` at κ=1; selected rows elsewhere; null rows of `repeated-groups.csv` with ≥5 rejections; design and calibration blocks of the three summaries |
| Fisher and allocation records | `score-design.csv` (all); all 36 event-count tables; `envelope-decomposition.csv` (all); selected projection and difference rows |
| Alignment reanalysis | `realized-envelopes.csv` (all); parts of `groups.csv` and `paired-differences.csv`; `gaussian-class-coverage.json` |
| Other | Known-pose entries of `paired-covariance-full-frequency-v1/summary.json`; three verification JSONs; `pose-metadata-inventory.json`; the focused view-variance audit and response |

**Not verified.**
- Array files were unavailable: score directions, calibration arrays, held-out scores.
- Execution: "25 targeted tests pass" and the replay counts are the authors' reports.
- Git ordering of declarations against outcomes.
- Third-party papers.
- Preferred-view and Monte Carlo raw JSON beyond their result reports.
- Older studies were rechecked only where Tables 1, 2, 4, 5, 6 meet their summaries. Those all replay.

## 1. Contribution and closest prior art

**Main text.**
1. A continuous L² audit of an affine estimator of a linear density feature, with sinc-kernel integrals, a quadrature remainder and a dual gap.
2. Proposition 1, with its sharp constant √3. It shows that a Gaussian prior with τ = B shares the fixed-pose uniform guarantee.
3. A conditional mixed-pose bound, which is vacuous in the one study that exercises it.
4. A post hoc diagnosis of alignment bias on 600 saved datasets, plus exact class-bias envelopes at realized poses.

**Appendices G–J.** A Monte Carlo goodness-of-fit test of a whole candidate map using power-spectrum and bispectrum contrasts. Its null allows amplitude in [0.9, 1.1] and a viewing density at most κ times Haar. It comes with several classical calibration bounds, a covariance-aware score and an allocation diagnostic.

**Closest supplied prior art.**
- Main text: Donoho (1994), Armstrong–Kolesár (2018), Backus–Gilbert (1968), Kanagawa et al. (2025), Fiedler et al. (2021), Ullrich et al. (2020), Rangan et al. (2024); for alignment bias, Jensen (2001), Shaikh et al. (2003), Henderson (2013).
- New branch: Zhang et al. (2024), which is the most direct since it compares candidate structures to image moments while profiling the viewing law; Sharon et al. (2020); Ortiz et al. (2020); Cossio–Hummer (2013) and CryoLike; Dufour (2006); Berger–Boos (1994); Rockafellar–Uryasev (2002); Duchi–Namkoong (2019); Maurer–Pontil (2009); Janon et al. (2014); Friedman (1989).

**Not in the packet (my own background, unverified here).** Kam (1980) on rotational invariants for randomly oriented particles; bispectrum inversion for multi-reference alignment (Bendory et al. 2018); Huber-type robust tests over density-bounded neighbourhoods.

**Net novelty.** Operator-specific integrals, one sharp elementary constant, and correct compositions of classical tools. No procedure outperforms a matched baseline in any tested regime. The authors say as much.

## 2. Technical checks

I re-derived each item and compared it with the code. All are correct as stated.

| Item | Note |
|---|---|
| Theorem 1; ridge stationarity λ = zb/(B²s); folded-width derivatives | q_s follows from homogeneity of q in (s, b) |
| Sinc Gram (7); cell transform | Factor ½ and sum-frequency term verified; cell sinc × d^(−3/2) matches `direct_cells` |
| Proposition 1, sharpness, adaptive-bias infimum | Reproduced the t⁴ coefficient φ(z)z(z²−3)/6 and 0.908 at c = 1.96 |
| Duality, feasible point, solve-gap identity (18) | |
| Quadrature remainder | Reproduced 3.2×10⁻²³ at radius 5, order 40 |
| Moments 1/12 and 13/720; remainder R₂; Theorem 2 | Valid; Theorem 2's premise is never met |
| Realized envelope (16) | Adjoint dephasing sign matches `dephase_observations` |
| Phase-model generating function and tail bound | Re-derived through the chi-square decomposition |
| G.1 variance identity | Needs a constant Laplacian; holds (4 per power feature, 0 for distinct-index triads) |
| G.3 binomial test; H.1 envelope, paired product, variance bound | The δ/4 + δ/4 + δ/2 split in the text matches `uq_view_variance.py:83-93` |
| H.2 CVaR identity, conditional Jensen, DKW stop-loss bound | |
| I.1 scale orthogonality; J constrained Rayleigh; J.2 convexity | Orthogonality holds on the training law only, as stated |

Seven findings change how Appendices G–J should be read.

**E. The test cannot detect its own design alternative on two of three candidates.**

The table uses the matched power score from the Fisher study: single-image exceedance frequencies at amplitude 1, and the null-envelope mean from the two allocations.

| Candidate | Null | 25% deletion | 50% | 100% | Envelope mean, L=32 / L=128 | Upper bound, κ=1 |
|---|---|---|---|---|---|---|
| 10028 | .4957 | .5104 | .5256 | .5578 | .5166 / .5137 | .5212 |
| 10049 | .5084 | .5242 | .5393 | .5658 | .5215 / .5185 | .5243 |
| 10076 | .5123 | .5163 | .5208 | .5283 | .5228 / .5198 | .5251 |

- On 10028 and 10076 the 25%-deletion frequency lies below even the point estimate of the envelope mean. The test rejects only when the frequency exceeds the bound, so its power tends to zero as particles are added. The same holds for all four score variants I checked.
- On 10049 the alternative sits at the bound. That explains projections of 4–9% at every sample size.
- Where the margin comes from, on 10028 (.4957 to .5212):
  - About 0.015 is the adversary that may pick amplitude within ±10% separately for each view. This is my estimate from the two allocations.
  - About 0.003 changes between L=32 and L=128.
  - About 0.005 is confidence slack.
- The 25% deletion moves the frequency by 0.015. The calibration-method and allocation studies (H.2, J.2) therefore address the two small terms.
- A plain z-test on the same contrast with no nuisance would need about 4,000 images on 10028 and 10049 for 80% power at 25% deletion, and about 50,000 on 10076 (squared separations .0013–.0017 and .00012 in `score-design.csv`). The nuisance envelope is the binding constraint, not the statistic's noise.

**F. Viewing-law tolerance is about 1–5% in χ², and nothing measures the real laws.**
- The worst-case shift of the null frequency is √(χ²·V), where V is the between-view variance of the conditional exceedance probability.
- For the scale-orthogonal power scores V ≈ .026, .0055 and .0014 (raw product estimates .0256, .0052, .0013).
- With 25% signals of .015, .016 and .004, the shift equals the signal at χ² ≈ 0.008, 0.05 and 0.014. A cap of κ = 1.1 already allows χ² up to 0.1.
- At κ = 2 every method is powerless on every candidate. On 10049 and 10076, full deletion at n = 10⁵ and the most favourable amplitude has projected rejection below 10⁻⁴⁰.
- No estimate of the stacks' viewing densities exists in the packet, although consensus poses for over 10⁵ particles per stack are on disk. The project's own note describes 10076 as having preferred orientations.

**G. Sensitivity to the noise level is unexamined.** This is a bound, not a verified number.
- Debiasing subtracts a known noise variance. A relative error ε shifts the score by ε·Σ_j w_j over the power block.
- The projection removes the mean-power profile, not the all-ones direction (`uq_candidate_score.py:36-38`).
- |Σw_j| can be as large as √220 ≈ 14.8, and the 25% design gaps are .012–.079. An error of a few tenths of a percent could equal the whole design gap. I could not compute Σw_j without the arrays.
- At least 80% of the Fisher separation comes from power features: adding the bispectrum raises it by 19%, 10% and 4%. With an unknown noise spectrum only the bispectrum block keeps a nuisance-free mean.

**H. The grouped envelope's extra assumption buys little.**
- Individual and grouped bounds differ by .0014–.0033 (for example .5195 against .5181 on 10028).
- The individual envelope also covers noise-dependent amplitudes. It could be used throughout at negligible cost.

**I. The authors' paired-variance bound is dominated where viewing variance is real.**
- On 10028 at κ = 1.1 the bounds are .5465 (split CVaR), .5495 (DKW), .5700 (paired variance) and .5699 (grouped ratio).
- The full-deletion projections there are .43 for split CVaR and 5×10⁻⁷ for paired variance.
- In the earlier oracle-score cases the paired bound was tighter only because its variance term was concentration slack, as the manuscript says.

**J. Low-order invariants discard most of the information (indicative).**
- The project's own known-pose record gives a mean squared mean-shift of .60, .96 and .16 per image at full removal. This is for the reference maps and their regions, at `probe_uq_full_covariance.py:55`. Scaled to 25% deletion it is .037, .060 and .010.
- The best moment contrast reaches .0016, .0017 and .00012. That is a factor of roughly 20–90.
- Caveat: the maps and regions differ between the two calculations. This is an order of magnitude, not a matched comparison.

**K. Null controls sit at the boundary, with one pattern to resolve.**
- The largest correct-null projection (.0537) occurs because the held-out null frequency (.58562) exceeds the calibrated bound (.58535). The cell is the 10049 raw power score at amplitude 0.9.
- In all six raw-score cells the held-out frequency at amplitude 0.9 exceeds the calibration estimate of the envelope mean, by .0002–.0032. In the population the envelope mean cannot be smaller.
- For the three independent power scores the differences are 1.4, 2.1 and 0.9 standard errors. This may be chance (about 2.5σ combined, post hoc). Only repeated calibrations can settle it.

## 3. Experimental audit

### 3.1 Kinds of evidence

| Level | What exists | What it supports |
|---|---|---|
| Conditional simulation coverage | 672 analytic checks; 18 true-pose fits replayed on 600 datasets | Implementation; representation bias is visible (true-pose noise-only coverage .92–.98) |
| Estimated-pose simulation | 200 independent datasets on each of three designs | Legitimate Monte Carlo trials per design, but for an aligner that degrades poses from 2° to 12–18° |
| Candidate goodness-of-fit in a simulator (new) | One training, one calibration and one held-out sample per candidate and study | Single-image frequencies conditional on one calibration; no unconditional size; no experimental images |
| Experimental prediction | Held-out image error; CryoLike scores | Unrelated to density coverage |
| End-to-end or experimental density coverage | None | — |

Three acquisition designs support no biological generalization. Results differ sharply across them. The asymptotic detection floor at κ = 1 lies between 25% and 50% deletion on 10028, between 10% and 25% on 10049, and near 50% on 10076.

### 3.2 Frozen versus development; independence

- Training, calibration and test use independent generator streams per stage and study (seeds 261016–261018 plus dataset plus stage×10⁶). I confirmed this in the runners.
- The allocation study inherits the Fisher scores and uses fresh calibration and test draws.
- Appendices G–J are sequential post-outcome development, as the authors state. Each study was designed after seeing the previous one.
- Scores, methods, amplitudes and deletions share draws. The 1,260 null cells are not independent trials.
- The candidates were fitted on the experimental stacks. A real application would need particles not used in fitting.

### 3.3 The candidate-validation branch: seven requested checks

1. **Whole candidate versus regional occupancy.**
   - The null is the entire candidate plus the imaging model. A rejection can come from any structural error, noise level, CTF, amplitude range or viewing law.
   - Held-out alternatives are exactly the designed deletion at the nominated region (`probe_uq_candidate_scores.py:79`). No off-region or global-error alternative is simulated, so specificity is unknown.
   - The region is the densest blob (20 Å SD, about 47 Å FWHM). Full deletion removes 9%, 27% and 15% of the candidate's L² norm.
   - The test is one-sided: it looks only for excess density in the candidate.
   - Any fitted map is wrong at the few-percent level, the same order as the detection floor. A point null with no equivalence margin is therefore ill-posed.
   - The manuscript's wording is suitably cautious. The experiments support no regional reading.
2. **Viewing law.** See finding F. The preferred-view controls supply κ exactly, use the oracle scores and are not the extremal tilt.
3. **Normalization and noise realism.**
   - One |CTF|/σ profile from the first particle is applied to every image (`probe_uq_bispectrum.py:48`).
   - Noise is white and known exactly. There is no defocus spread and no window-induced correlation.
   - The noise scale is the supplied scalar. Earlier rounds found the experimental calibration could bound the noise SD only within 1.9–4.6 times that value. I did not recheck this.
4. **Training, calibration, test independence.** Sound within the simulator (3.2). Unconditional size is untested because each study has one calibration per candidate.
5. **Max-over-noise versus noise-adaptive amplitude.** A second-order effect here (finding H). The first-order term is view-dependent amplitude (finding E). Whether real amplitudes vary with view by ±10% is unmeasured.
6. **Comparator selection.**
   - The classical bounds are appropriate and retained, and no unadjusted minimum is taken.
   - Missing: any information ceiling, a likelihood-ratio baseline such as the CryoLike machinery the authors already ran, and the moment metric of Zhang et al.
   - Figure 8 shows only the authors' procedure.
7. **Power-interval interpretation.**
   - "Projected rejection probabilities" are binomial transforms of one held-out frequency given one calibration.
   - At n = 10⁵ the transform is so steep that a band spans [.16, .99] (10028, 50% deletion).
   - A .006 change in the bound moves a projection 45-fold (10076: .000176 to .0080).
   - The Fisher-minus-matched intervals are Bonferroni-wide ([−.59, .74]), though the retained joint tables would allow a paired analysis.
   - The only observed test outcomes are the 128 groups of 1,024 images.

### 3.4 Main text since round 3

- The bias finding is now reported, with Figure 1 and Table 6. Its numbers replay: .720 against .955; shifts of .357 ± .029 and .325 ± .024; projected noise of 1.23 and 1.31 SDs.
- Realized envelopes replay (.306–.532 of no-data). They are computed only at the realized 12–18° errors, so they say nothing about realistic radii.
- The direct folded-width search gains .06–.10%. This confirms the audit is a point on the ridge path.
- No redesigned simulation with a working aligner, no latent-pose test of Theorem 2, no lower bound and no new construction.

### 3.5 Usefulness and scope

- On 10049 the audit estimator's RMSE is 1.39–1.60 times the pilot's own error for the central feature, even at true poses.
- No interval excludes zero under any rotation budget.
- The candidate test's best case is detecting removal of a quarter of a molecule's L² norm, under exactly uniform views.

### 3.6 Statements not supported as written

1. **Abstract, "re-estimates local poses".** The headline .72 coverage belongs to a refinement that worsens poses sixfold. The abstract does not say so.
2. **§6.4, "even the realized oracle envelopes are large".** True only at 12–18° errors.
3. **Appendix I, "reveal sensitivity limits".** The evidence is stronger: the design alternative lies inside the null envelope (finding E).
4. **Appendix J.2.** The allocation question concerns about 0.003 of a 0.021 margin.
5. **"Candidate-derived".** The score needs no reference map, but every evaluated alternative is the designed one.

## 4. Concerns

### Status of earlier concerns

| ID | Concern | Status |
|---|---|---|
| R1 | Experimental inputs assumed | Open |
| R2 | No useful inference | Open; see 3.5 |
| R7 | Novelty of the audit | Open; ridge equivalence now confirmed by the authors' own search |
| R9 | Theorem 2's premise never met | Open |
| R11 | No informative end-to-end study | Open; no working-aligner study run |
| R13 | Bias finding unreported | Resolved |
| R14 | Pose bound structurally vacuous | Partly: diagnostic done; no theorem or new construction |
| R15 | Density class and metric | Open |

### Major concerns

**R16 (new). The manuscript is two unrelated half-papers.**
- *Evidence.* `focused-main.tex` contains no reference to Appendices G–J. The estimand changes from an interval that is uniform over a density class to a point-null test that assumes the fitted candidate exactly right.
- *Change.* One paper, one estimand. Remove G–J or rebuild around them.
- *Acceptance.* Abstract, introduction and discussion cover every method evaluated.
- *Status.* Fixable by restructuring. Neither half is acceptance-worthy alone.

**R17 (new). The candidate test is structurally insensitive.**
- *Evidence.* Finding E.
- *Change.* Produce the ledger in step 1 of the plan. Then measure the amplitude nuisance or change the statistic.
- *Acceptance.* The 25% alternative lies outside the population null envelope on at least two candidates, under nuisance sizes measured on the stacks.
- *Status.* Fatal to this construction. No score, bound or allocation variant fixes it.

**R18 (new). The viewing-law model is unrealistic and unmeasured.**
- *Evidence.* Finding F.
- *Change.* Measure the stacks' viewing densities. Replace the Haar cap by an estimated law with quantified error, or use per-particle pose information.
- *Acceptance.* The measured law, with split-half stability, lies within the tolerance at which the 25% alternative remains detectable on at least two stacks.
- *Status.* Fatal as posed.

**R19 (new). Simulator realism; no experimental contact.**
- *Evidence.* 3.3 item 3; finding G.
- *Change.* Per-particle CTFs, an estimated noise spectrum with its uncertainty, and real-data controls on the held-out exposure cohorts.
- *Acceptance.* On real held-out particles, the intact fitted candidate is not rejected and a fully deleted candidate is rejected, on at least two stacks.
- *Status.* Essential if the branch continues. Outcome uncertain.

**R20 (new). Interpretation and specificity.**
- *Evidence.* 3.3 item 1.
- *Change.* Drop all regional language and present an omnibus compatibility test with an equivalence margin, or add off-region and global-error alternatives.
- *Acceptance.* Rejection at most .1 for an equal-L² perturbation elsewhere, or no regional claim.
- *Status.* Fixable.

**R21 (new). No unconditional validity evidence.**
- *Evidence.* Finding K; 3.3 item 7.
- *Change.* At least 200 independent repetitions of calibration plus test, with actual datasets, at the boundary null and under an extremal tilted law.
- *Acceptance.* Realized size at most .05, with a 95% upper limit at most .08.
- *Status.* Fixable. Worth running only if R17–R19 are met.

**R22 (new). Comparators and novelty of the branch.**
- *Evidence.* Findings I and J. Every ingredient is classical, as the text says.
- *Change.* Report the known-pose and Haar-marginal likelihood-ratio ceilings and one likelihood baseline on the same simulated images.
- *Acceptance.* The procedure is within a factor of three of the marginal-likelihood ceiling in squared separation, or beats a likelihood baseline under a stated nuisance.
- *Status.* Not fixable by wording.

**Carried major concerns.**

| ID | Needed | Acceptance | Status |
|---|---|---|---|
| R11 | One simulation with an aligner that improves on initialization; true-pose, same-image and independent-image arms; one in-class adversarial generator | Gate passes; the bias finding persists or is withdrawn | Essential if alignment claims are kept |
| R14 | A lower bound, or exact envelopes at 0.5–2° and a construction built on them | A theorem, or widths below 0.5 of no-data and within twice fixed-pose | Fatal to the pose-aware contribution as formulated |
| R9 | A latent-pose model in which Theorem 2's premise holds by construction | Nominal coverage in the R11 study | Potentially fixable |
| R7, R2, R15, R1 | As in round 3 | As in round 3 | Open |

### Minor concerns

- **M1 (carried).** Not anonymized: the repository URL and a self-citation identify the author.
- **M8 (carried).** Low (1997) is uncited; Cai–Low is in the bibliography file but not cited.
- **M25 (carried).** The mixed interval does not reduce to Theorem 1 as radii vanish (multiplier 2.72 against 1.96).
- **M28 (carried).**
  - Page 5 has a blank band under Table 2.
  - Pages 19, 20, 22 and 23 are largely white.
  - Three of six panels in Figure 9 are flat at zero.
- **M31.** The paper contains review-process commentary ("three full-paper rejections remain unresolved", "no acceptance reassessment has occurred"). A paper is not a lab log.
- **M32.** Appendix J.2 gives .000182 for the 10076 combined-Fisher split-CVaR projection (`fisher-score-development.tex:129`). The result CSV and the project's own table give .0001755.
- **M33.** Explain the .0537 null projection and the pattern in finding K.
- **M34.** Figure 8 should show the best classical comparator beside the paired-variance curves.
- **M35.** Cite the invariant-feature and robust-testing literature noted in section 1.
- **M36.** `CURRENT-EVIDENCE.md` links a round-3 response plan that is not in the evidence directory.
- **M37.** Appendix F and E.7–E.8 are tangential to either half.

## 5. Revision plan

Three constructions in a row have failed their own gates: worst-case pose balls, pose-invariant paired statistics, and bounded-view moment tests. The common cause is treating poses, views and amplitudes as adversarial nuisances to be bounded or removed by invariance. That discards the per-particle information that makes cryo-EM work, and what remains is swamped by the nuisance class. The direction change is to model these quantities as latent, measure their distributions on the data, and calibrate statistically with real-data controls.

**Essential, in order.**

1. **Freeze variants and build a ledger (days; existing arrays).** No new score, bound, allocation, κ or sample size. For each candidate and deletion level report:
   - the per-image squared separation for the known-pose likelihood ratio, the Haar-marginal likelihood ratio and the moment contrast;
   - the null-envelope decomposition (fixed amplitude, unknown global amplitude, view-dependent amplitude, noise-adaptive amplitude, finite-L bias, confidence slack), with the asymptotic detection floor at each stage;
   - the viewing-law and noise-level tolerances; Σw_j; a bispectrum-only separation.

   This resolves whether weak power comes from the data, the compression or the nuisance model.
2. **Measure the nuisances on the three stacks (days; existing metadata).** Viewing-density ratio and χ² against uniform, with split-half stability. Amplitude dispersion and its dependence on view. Noise-spectrum variability. Defocus spread. This resolves whether κ ≤ 1.1, ±10% amplitude and exactly known white noise are within an order of magnitude of reality.
3. **Decide once.**
   - If the marginal-likelihood ceiling is at least ten times the moment separation, as finding J suggests, abandon moment contrasts. Test regional occupancy with a latent-pose likelihood statistic, nuisances estimated per particle, and Monte Carlo calibration. Per `CURRENT-EVIDENCE.md`, which I did not re-read in detail, an earlier certified-mixture attempt left global gaps above 12,000. Certificates are therefore not the route.
   - If even the ceiling cannot give 80% power at 25% deletion with 10⁴ particles under the measured nuisances, the question is unanswerable at this scale. Write that up as the result.
4. **One pre-registered gate for whatever survives.** All of:
   - size as in R21;
   - power at least .8 at 25% deletion on two of three candidates under measured nuisances;
   - specificity as in R20;
   - the real-data control in R19.
5. **For the main-text half, if kept.** The working-aligner study (R11), exact envelopes at realistic radii, and a lower-bound theorem (R14).
6. **Rewrite as one paper (R16).**

**Stopping rule.**
- A variant may run only if it targets the largest term in the step-1 ledger and comes with a written prediction that it moves the detection floor below the gate.
- Two consecutive variants that miss their prediction end the branch.
- If the gate fails after one revision, end the branch.

**Optional.** Experimental density intervals with evidence-based inputs (R1); CTF and scale nuisances in the audit; additional stacks.

**Not needed.** More calibration bounds, allocations, score directions, cap laws, κ values or particle counts. More radius-5 replicates of the current aligner. Survey growth.

## 6. Verdict

**Reject. Confidence 4 of 5.** I read the theory, the code paths behind the new headline numbers and a broad sample of raw records. I ran nothing and could not open arrays.

**Is the current work a strong contender for ICML acceptance? No.** Reasons, in order of weight:

1. No proposed procedure is informative in any regime the evidence supports. The pose-aware intervals are vacuous, and the candidate test cannot detect its own design alternative on two of three candidates.
2. The fixed-pose audit is matched by a Gaussian baseline with the same guarantee. The new branch's one original bound is dominated by classical CVaR where it matters.
3. The new branch assumes near-uniform views, exactly known white noise and a single CTF. None of this is measured, and no experimental image is tested.
4. The main empirical finding comes from a degenerate aligner. The redesigned study requested in round 3 was not run.
5. The manuscript is two unconnected halves, with a development log as appendix.

**Worth preserving.** The sharp envelope proposition, the matched baseline, the honest null and failure reporting, the exact realized-envelope machinery, and the replay discipline. More experiments, correct classical algebra and open records are real virtues, but they are not a contribution.

**What would change the assessment.** A procedure that passes the gate in step 4 on real held-out particles. Alternatively, a theorem showing why no honest procedure can, paired with a clean estimated-pose study.
