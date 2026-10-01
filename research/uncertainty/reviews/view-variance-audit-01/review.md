# Focused audit: bounded-view candidate validation (paired-power-v1, view-variance refinement)

**Bottom line.** I found no error that invalidates the stated conditional guarantee in exact arithmetic. I found one unstated assumption (V1), one missing sharp comparator that bears directly on the contribution (V3), and back-calculated evidence that the headline view-variance numbers are set by concentration slack rather than by measured viewing variance (V4). The negative conclusions stand: no useful power at kappa ≥ 2 anywhere, and none on 10028 even at kappa = 1.

## Findings

### V1 — Major (assumption gap): the grouped envelope needs amplitude ⟂ noise given view

**Location:** VIEW-VARIANCE-PROTOCOL.md §"Amplitude envelopes", sentence "This holds even if amplitudes vary across actual particles and depend on orientation"; docstring of `grouped_amplitude_envelope`.

**What checks out.** For fixed `a` in cell `j(a)`, `1{f(a m+eps_l)>c} <= Z_{l,j(a)}` for every `l`, so the replica average is `<= X` pointwise. Taking the noise expectation gives `eta(R) >= sup_a Pr(f(a m(R)+eps)>c | R)`. A random amplitude with a law depending on `R` but independent of the particle's noise is then covered by integrating over `a`.

**What is missing.** Domination holds only for an amplitude that does not select on the realized noise.
- The per-noise envelope `T` (max inside each noise draw) covers even noise-adaptive amplitudes.
- The grouped envelope `X` does not: an amplitude correlated with noise can have event probability up to `E[T|R] > eta(R)`. Examples are ice thickness driving both, picker or 2D-class selection, and per-particle scale estimated from the image.
- So `grouped_ratio` and `view_variance` protect a strictly smaller null than `individual_ratio`. Part of the individual→grouped gain is bought by assumption, not by statistics.

**Fix.** State "amplitude conditionally independent of noise given view, and independent across particles" in the theorem. Label the three methods by the null each protects. Falsification test: simulate `a_i = argmax_a f(a m+eps_i)` and confirm the grouped bound is violated while the individual bound holds.

### V2 — OK (verified algebra): centered product, variance subtraction, minimum, binomial domination

**Paired independence.** The script draws noise with shape `(n,2,REPLICAS,nq)` in one iid call and both groups share the same `mean`. `X1` and `X2` are therefore conditionally iid given `R`.

**Center.** `b` comes from the older run. Stage-0 seeds were 271040 and 271061 there, against 271041 and 271062 here, so the streams are distinct. With `b` fixed, `E[(X1-b)(X2-b)|R] = (eta_L(R)-b)^2`.

**Variance subtraction.** `E_Q Y = V + (mu-b)^2`. On `{mu in [L_mu,U_mu]}`, `(mu-b)^2 >= dist(b,[L_mu,U_mu])^2`, so `V <= U_Y - dist^2`. This is the correct one-sided direction. The `min(1/4)` and `max(0)` clips are valid.

**Range of Y.** `[-b(1-b), max(b^2,(1-b)^2)]` is exact for `X in [0,1]`.

**Empirical Bernstein.** Maurer–Pontil Theorem 4 is one-sided with `log(2/delta)`, a `7/(3(n-1))` term and the unbiased variance, as I recall it. The code matches that form, and the affine rescaling is right. The allocation `delta/4 + delta/4 + delta/2 = delta` matches between protocol and code, and it allows arbitrary dependence between the three events.

**Cauchy–Schwarz.** `E_P eta = mu + E_Q[(w-1)(eta-mu)]` and `E_Q(w-1)^2 = E_Q w^2 - 1 <= kappa-1`. This is correct given `P << Q`.

**Minimum.** All three expressions are deterministic upper bounds on `E_P eta`, monotone in `(mu,V)`, so the minimum holds on the single joint event without correction.

**Binomial domination.** Conditional on calibration, the indicators are independent with `p_i <= p_bound`, and uniform coupling gives domination by `Bin(n,p_bound)`. Heterogeneous `P_i` and amplitude laws are fine. The `rejection_count` search is monotone and returns `n+1` at `p=1`. The level `.049 = alpha - delta` is right.

**Two caveats to write down.**
- `eta` is `eta_L` with L = 32, so `V` is the variance of the max-biased envelope mean, not of the true sup-probability. That is valid but should be explicit.
- Independence of indicators is on the tested set. Any image-dependent particle selection breaks both the noise law and V1.

### V3 — Major (design and novelty): the Cauchy–Schwarz bound relaxes a known exact answer, and the comparator is missing

**Location:** PROTOCOL "E_P eta <= mu + sqrt((kappa-1)V)"; "Compare three prespecified methods"; PRIOR-ART.md.

**Reasoning.** Under the stated assumption `0 <= w <= kappa`, `E_Q w = 1`, the exact worst case is classical:

    sup_P E_P eta = CVaR_Q at tail mass 1/kappa of eta   (mean of the top 1/kappa fraction)

It is attained at `w = kappa·1{eta in top 1/kappa}`. Three consequences:
- It is never worse than `kappa·mu` or the Cauchy–Schwarz bound at the population level.
- By Jensen, `CVaR(eta) <= CVaR(Xbar)` with `Xbar = (X1+X2)/2`. A valid finite-sample bound therefore follows from the empirical tail mean of the 32,768 view means plus a DKW or CVaR concentration term. No paired product and no center `b` are needed.
- The Cauchy–Schwarz bound uses `kappa` only through `chi^2(P||Q) <= kappa-1`. The honest assumption for that method is a chi-squared budget, which is weaker and more estimable than a sup-norm ratio.

Without the CVaR comparator, "less crude than `kappa·mean`" is a comparison against the weakest baseline. Two cheap ablations are also absent: `V <= U_mu(1-L_mu)`, and `Var(Xbar) >= V` without pairing.

**Fix.** Add the CVaR-of-`Xbar` bound and both ablations on the existing arrays, under a newly frozen protocol. If CVaR matches or beats `view_variance`, the paired-product machinery is not the contribution. Restate the variance method under a chi-squared budget.

### V4 — Major (interpretation): view-variance results look slack-dominated

**Location:** RESULTS table, `view_variance` rows.

**Reasoning.** This is my own normal-approximation back-calculation from the table, not from the arrays.
- On 10049/power, the kappa = 1 → 1.1 collapse of `grouped_ratio` (0.561 → 2e-17) implies `mu ≈ 0.4`.
- The `view_variance` values 0.546 / 0.116 / 3.2e-5 / 3.4e-16 are then mutually consistent with a single `V_U ≈ 4e-4`.
- At `M = 32,768`, `delta/2 = 5e-4` and `b ≈ 0.4`, the Bernstein range term alone is `7(u-l)log(2/d)/(3(M-1)) ≈ 3–4e-4`, and the sampling term adds about 1–2e-4.

So `V_U` is roughly all concentration slack, and the true `V` under Haar is plausibly near 0. Two implications:
- The experiment has not exercised a case where viewing variance matters. It shows only that these scores have nearly view-independent null event rates under Haar.
- The kappa ≥ 1.1 numbers are governed by `sqrt((kappa-1)·O(1/M))`, a budget artifact rather than a property of the specimen.

**Fix.** Report `centered_product_mean`, `centered_product_upper` and `joint_mean_interval` (already saved) with a decomposition of `V_U` into estimate, sqrt term and range term. If the slack is confirmed, say so, and do not describe the result as a measured viewing-variance sensitivity.

### V5 — Moderate (fairness): what "matched budget" means

**Location:** PROTOCOL "same simulation budget"; `viewing_probability_bounds`.

- All three methods share 32,768 projections × 64 noises. This is fair in projection cost.
- `individual_ratio` is not the v0.7.3 procedure (131,072 views, Clopper–Pearson). It is an empirical-Bernstein bound on view means with an additive `≈5e-4` range penalty, at a split chosen for the variance method.
- Ordering between methods is deterministic given `p_bound`, because all share one held-out count. The overlapping pointwise intervals therefore neither support nor undermine the ordering. The missing uncertainty is calibration-replicate variability, since there is one calibration draw per stack.
- The grouped gain at kappa = 1 is small on 10049 (0.505 → 0.561) and partly attributable to V1.

**Fix.** Add a best-allocation baseline for the ratio methods at equal projection count. Bootstrap or re-seed calibration to show the spread of `p_bound`. Present paired differences of critical values rather than overlapping intervals.

### V6 — Moderate (conclusions supported): narrow and largely negative

**Location:** RESULTS "Interpretation".

- **Supported:** `X <= T` pointwise. The variance bound yields a smaller `p_bound` than `kappa·U` at kappa = 1.1 in all four cases.
- **Not supported as usefulness:**
  - At n = 10,000 and kappa = 1.1, the best power is 0.185 [0.075, 0.36].
  - 10028 is ≤ 0.008 at kappa = 1.1 and only 0.02–0.19 at kappa = 1. `individual_ratio` on 10028/power is below alpha.
  - Power at kappa ≥ 2 is negligible everywhere.
- **Conditions:** oracle direction and threshold, full-region removal, Haar alternative, amplitude exactly 1 in the alternative, one transfer profile, four post-outcome selected cases.
- **Presentation:** entries printed "0 [0, 0]" are floating underflow and should read "<1e-300". The true-map type-I projections are not in the packet's table and should be shown beside power.
- **Alternative realism:** the held-out alternative fixes amplitude at 1 while the null is enveloped over `[.9,1.1]`. Power under amplitudes spread across the interval is unreported and could be lower.

### V7 — Minor (numerics): cubic and cell edge cases look sound, uncertified as disclosed

**Location:** `cubic_cell_events`, `cubic_interval_maximum`.

- **Checked by reading:**
  - The derivative is `a x^2 + b x + d` with `a=3c3, b=2c2, d=c1`, solved by the stable `z` form.
  - `z == 0` occurs only for the double root at 0, which `z/a` handles.
  - The power-only case has `c3 == 0` exactly and takes the linear branch.
  - Marking a cell true for any in-cell stationary value above threshold is correct, including minima.
- **Residual risks:**
  - A discriminant rounded slightly negative drops a near-double root. That is an inflection, so it is harmless.
  - A stationary point within an ulp of a cell edge may be assigned to the neighbouring cell. The value error is second order, but it is not certified.
  - A `-0.0` in `b` flips `copysign`, which is harmless.
  - Tiny nonzero `a` can produce a huge or infinite root, which the validity test rejects.
- **Test gap:** tests use generic normal coefficients. There is no near-degenerate case, and no coverage test of the full bound.

**Fix.** Add adversarial cubic tests: root exactly at an edge, `disc ≈ 0`, `|c3| ~ 1e-300`. Add a Monte Carlo coverage test of `viewing_probability_bounds` with a known `eta` and the extremal `w`.

### V8 — Moderate (verification scope): the replay is narrower than "independent"

**Location:** `verify_uq_view_variance.py`; provenance JSON.

- The reported counts are internally consistent: 8×64 = 512 polynomials, 2·2·2·5·3 = 120 bounds, 360 critical values, 8 held-out counts.
- The replay re-derives polynomials by Vandermonde interpolation through `moment_features` for one view per stack and stage. The other 32,767 views rely on the production path.
- It reuses the stored `example_direct` mean, the production `moment_features`, and the same Bernstein formula. It therefore checks transcription, not validity.
- I could not confirm that `data['transfer']` is the "first saved transfer profile", because the script uses the array as-is.
- I could not confirm whether the frequency set `q` contains antipodal pairs. If it does, Hermitian symmetry of real images contradicts the independent complex noise per coordinate.

**Fix.** Replay a random ~1% of views drawn with a verifier-chosen seed. Assert `transfer.ndim == 1` and record its provenance. Assert that no `q_i = -q_j`, or model the Hermitian coupling.

### V9 — Major (relevance): the sup-norm ratio to Haar is the wrong scale for cryo-EM

**Location:** both protocols, assumption `dP_i/dQ <= kappa` with `Q` = Haar.

Real orientation distributions are concentrated, so I would expect the sup-density ratio to Haar to be far above 5 in many datasets and not estimable without smoothness assumptions. The packet contains no measurement of this. All methods have zero power by kappa = 2–5. As posed, the guarantee is vacuous exactly where preferred orientation exists.

**Fix.** Take `Q` to be an estimated viewing law, held out or from a consensus refinement, and assume a small chi-squared misfit (see V3). Calibrating that misfit is then the open problem.

## Could this be a meaningful contribution, and what is classical

**Classical ingredients:** the Monte Carlo test with nuisance envelope, Clopper–Pearson and empirical Bernstein bounds, binomial domination, the pick-freeze covariance identity, chi-squared and Cauchy–Schwarz robustness, and the CVaR representation of the bounded-ratio worst case. The composition is correct but not a theorem-level contribution, as the packet itself says.

**What could be a contribution** is an empirical, problem-specific fact, if it holds: that for translation-invariant moment scores the null exceedance probability is nearly view-independent (V4). That would make view-law robustness cheap. It is currently an inference from slack, not a demonstrated result, and it sits behind an oracle score and assumed noise, CTF and amplitude.

## Most discriminating next experiment

Run the adversarial extremal tilt, not a generic preferred-view model.

1. For each kappa, draw null data from `w* = kappa·1{eta in top 1/kappa}`, with `eta` estimated on an independent split.
2. Run real repeated datasets, not binomial projections, and report realized type-I against each method's `p_bound`.
3. Report the CVaR comparator and the `V_U` decomposition alongside.
4. Add a second arm with the alternative also tilted, and a score or region chosen so that `eta` varies materially with view.

Possible outcomes:
- CVaR ≈ `view_variance` and realized type-I far below both: the method is slack-limited, and the contribution reduces to "eta is flat".
- Any realized type-I above alpha: a bug or a violated assumption (see V1).
- Power vanishing once `V` is non-negligible: the refinement does not help where it would be needed.

## What I did and did not verify

**Verified by reading and hand derivation only:**
- The domination, product-moment, variance-subtraction, Cauchy–Schwarz, joint-event and binomial-coupling arguments.
- Agreement of the code with the protocol on the delta split, ranges, the level `.049`, seeds, batch sizes and array shapes.
- The quadratic-root branch logic.
- Arithmetic consistency of the 512/120/360/8 counts and the table's monotonic patterns.
- A rough normal-approximation back-calculation of `mu` and `V_U`. It is approximate and must be confirmed from `summary.json`.

**Not verified:**
- I executed no code and inspected no arrays, CSV or JSON beyond the pasted text.
- I did not read Maurer–Pontil, Duchi–Namkoong, Janon et al., Dufour or Berger–Boos. The Theorem 4 form and the CVaR duality are from memory.
- I did not see `moment_features`, `CellObservationOperator`, `direct_cells`, the transfer shape, the `q` geometry or the old Monte Carlo summary.
- I did not see the true-map type-I projections, or any n = 1,000 / 100,000 or kappa = 1.01 results.
- Nothing here bears on experimental calibration, novelty, or the three standing rejections.
