# Post-round-4 methods consultation

**Basis.** I read only the supplied text. I ran nothing, opened no arrays and read no outside sources. Figures marked *(est.)* are my hand calculations from the summary tables and need replaying from the saved arrays. Literature cited from memory is marked *(not in packet)*.

## Summary

- **Diagnostics.** They close the moment branch and show that a likelihood-type statistic has 29–55 times more separation. They never expose that statistic to a nuisance, so they do not yet justify a method.
- **Estimand.** One candidate merits a bounded test: a population fraction, treated as a linear functional of an unrestricted joint state–view law. Component occupancy is the special case where state B is state A without the region. The evidence justifies a three-day arithmetic gate (G1), not the method build, and G1 may well fail for the current compact regions.
- **Integration.** The catalogue revision is sensible repair but off the critical path. Run it once, after a free autopsy of the saved draws, then close the numerical branch whatever the outcome.

## 1. Do the diagnostics answer review 4?

| Review-4 question | Status | Evidence |
|---|---|---|
| Is weak power due to data, compression or nuisance? | Answered for moments | Fixed-bank likelihood-ratio D is 29×, 40×, 55× the best moment D. At 25% deletion that is about 130, 90, 920 images for 80% power, against 4,000–50,000 |
| Is κ ≤ 1.1 near reality? | Answered: no | Recorded cross-half χ² of 1.07, 2.0, 0.25 against tolerances of 0.008, 0.05, 0.014 |
| Is known white noise tenable? | Partly | Moments need 0.26–0.5% accuracy on two stacks; corners show a 3.8–5.7-fold power range; the in-band spectrum is still unknown |
| Does a likelihood statistic survive measured nuisances? | Not addressed | D3 |
| Unconditional size; experimental control | Not addressed | Correctly deferred |

### Diagnostic issues

**D1. The tenfold trigger is met without converged integration.**
- A log ratio computed on a fixed orientation bank is itself a statistic, so its measured D is an achieved separation and the best statistic can only do better.
- The two banks agree in D to under 1%, and null variance is about 1.09·D on all three stacks, as a likelihood ratio should give.
- The protocol's sentence that an unstable approximation cannot "trigger a claimed tenfold information advantage" is therefore too cautious in this one direction.
- Non-convergence does block calling the number a KL value or ceiling, any impossibility claim, and any use that compares likelihoods across parameters.
- Quote the sampling standard errors from the JSON beside D.

**D2. On 10028 the "Haar-marginal" column is a best-template statistic.**
- A median ESS of 1.23 means one bank point dominates, so label it as a fixed-bank template ratio.
- The losses against known pose (7%, 27%, 62%) mix real pose uncertainty with grid loss, and the split is not identified.

**D3. The ledger is asymmetric; this is the largest remaining confound.**
- Moments were stressed against amplitude and noise. The likelihood statistic was evaluated only at Haar views, one CTF, amplitude one and known white noise.
- Missing rows, all computable from existing arrays:
  - **Amplitude-orthogonal information**, `δ²·E[|g|² − |<g,m>|²/|m|²]`. A deleted region and a lower amplitude differ only through the part of the region's projection not collinear with the particle's.
  - **Null-mean shift** of the statistic under the measured non-flat background shape and under a recorded non-Haar law.
- With known pose the linear log ratio has null mean `−|Δ|²/2` for any zero-mean noise, so noise enters only through pose inference. That suggests far lower sensitivity than the moments, but it is unmeasured.

**D4. The event decomposition is sufficient.**
- At κ = 1.1 the view-cap increment alone is 3.5, 1.5 and 2.8 times the signal.
- No nested calculation of inner-replication bias is worth running.

**D5. Both nuisance classes in the frozen test were mis-sized, in opposite directions.**
- The view cap was too small by one to two orders of magnitude in χ².
- The ±10% view-dependent amplitude looks too generous: the recorded alpha–view association is 2.6% on 10028 and 0.74% on filtered 10076. Alpha's meaning is uncalibrated, so this is indicative only.
- The 10076 filter effect is the only supplied evidence that a viewing law can differ between subpopulations.

**D6. Corner patches cannot deliver the in-band noise spectrum.**
- An 8×8 patch resolves frequency in steps of one-eighth of the sampling rate, with the mean removed and particle tails present.
- The retained low frequencies are then constrained by a handful of leaky coefficients. I cannot check the band geometry from the packet.
- R19 needs noise-only windows.

### Numerical issues

**N1. The failure signatures differ by stack *(est.)*.**
- **10028:** the ESS fraction of 0.552 equals the ceiling I derive for one Gaussian mode under this proposal.
  - In tangent coordinates the ACG is a 3-D Cauchy with the Laplace scale, giving `∫p²/q = (√π/8)(1+3+3.75) = 1.717`.
  - With the 95% share the ESS fraction is 0.553, so the eight-start machinery contributes nothing there.
- **10049:** ESS is fine, but the between-bank RMS is about five times what the median ESS predicts. The failure sits in a minority of images, which points to missed or mis-weighted modes.
- **10076:** fails both criteria, with RMS about twice the median-ESS prediction.
- These use a posterior SD of the log ratio of about 0.13 (10049) and 0.09 (10076), inferred from the uniform-bank results.

**N2. The two gate criteria are not equivalent *(est.)*.** Having 90% of images within 0.01 needs a lower-decile ESS near 900 on 10049 and 450 on 10076. A median of 256 is much weaker.

**N3. Ratio agreement leans on common-mode cancellation.**
- The RMS log-integral discrepancy (0.015, 0.15, 0.16) is 7–20 times the ratio discrepancy.
- 10028 passes for the ratio only.
- Any method needing weights across view components or distinct states inherits the larger number.

**N4. The conditional-unbiasedness note is correct; three sharpenings.**
- **Variance.** With a bounded integrand and 5% Haar mass the variance is always finite; its size is the issue. A 2° mode has Haar mass near 10⁻⁵, so the roughly 410 defensive draws per bank find a missed mode with probability under 1%.
- **Products.** Independent unbiased per-image estimates give an unbiased stack likelihood.
  - A split-likelihood-ratio numerator therefore needs no bound at all.
  - One Markov step gives `Σ log I_i ≥ Σ log Î_i − log(1/δ)`.
  - Only the null upper bound needs a certificate, which is the hard half of the note's two-sided requirement.
- **Jensen bias.** About `−1/(2·ESS)` per log integral, or 0.002 at ESS 230. It matters when ESS differs between hypotheses.

**N5. The rare-mode example is correct, and simulation allows a stronger necessary check.**
- The integrals are `p/(1+p)` and `2p/(1+p)`, and the probability bound holds.
- Here the proposal is fitted to the midpoint map and both maps share modes, so a missed mode largely cancels in the ratio. The exposure is in absolute integrals (N3).
- **Planted truth.** In simulation the generating pose is an exact posterior draw under the matched model.
  - Replace one draw by it, and `1/Î_planted` is unbiased for `1/I` (the bidirectional Monte Carlo identity, *not in packet*).
  - That gives a per-image missed-mass detector and a stack-level bound `Σ log I_i ≤ Σ log Î_planted,i + log(1/δ)`.
  - The protocol permits generating poses for reporting. No such check exists for experimental images.

**N6. Necessary gate versus sufficient certificate.**
- **The gate shows:** Monte Carlo noise of one ratio, at two maps, in the simulator, given a proposal.
- **It cannot show:**
  - common misses under a shared proposal;
  - uniformity over a parameter set;
  - a joint budget, where systematic per-image error b must satisfy `b ≪ sqrt(I/n)` for estimation and `n·b ≪ 1` for product-likelihood e-values;
  - posterior geometry under experimental model error.
- **A sufficient certificate** needs joint upper bounds on unseen mass. Earlier certified attempts left gaps above 12,000, and I do not recommend another.
- The better route is an inference in which integration error costs efficiency, not validity (E3).

**N7 (minor).**
- Only the first of 32,768 regenerated rotations is asserted against the saved training means. Assert a few later-batch indices; the 10028 result suggests the replay is right.
- Mode weights use Laplace mass with clipped covariances. Check the failing images for clipped eigenvalues.

## 2. Estimand and method

### Choice of target

| Target | Verdict | Reason |
|---|---|---|
| Regional density confidence | No | Three constructions failed; nothing survives latent poses; no truth exists on the stacks |
| Omnibus compatibility | Not as headline | With 10⁵ particles and this statistic every fitted map is rejected; an equivalence margin needs an effect size, which leads to the next row |
| Population functional | Only candidate | Linear in the mixing law, so convex; yields an interval; the viewing law is profiled exactly, not capped |

**E1. Estimand.**
- Each particle has latent `ξ = (z, R, a)` drawn from an unknown law μ on `{A,B} × SO(3) × [a₋, a₊]`.
- μ is unrestricted: no Haar assumption and no state–view independence.
- The whitened image is `y = a·C_i·P_R·V_z + noise`, and the target is `π = μ(z = A)`.

**E2. What pooled reweighting estimates (falsifiable).**
- With known poses the conditional score mean is `(π_R − π)·J(R)`, where `π_R` is the local state fraction and `J(R) = ∫(p_A − p_B)²/p_mix`.
- The pooled estimate therefore converges to `π*`, with `π*/(1−π*) = π/(1−π) × E_νA[J]/E_νB[J]`.
- Pooled inference is biased exactly when the two states' viewing laws differ in mean discriminability.
- In the weak-signal limit J is about `|g(R)|²`, which the saved region means give directly.
- This is textbook for known poses. The latent-pose version is open.

**E3. Method: a view-unbiased calibrated score.**
1. Fix a bounded feature map `h(y)`, for example soft template-bank assignment to state × view cell. Any numerical approximation is part of h.
2. Estimate the response `k(ξ) = E[h(Y) | ξ]` by simulation with the same code.
3. On a separate fold choose λ to minimise variance subject to `|λ₀ + λ·k(ξ) − 1{z=A}| ≤ β₀` at every catalogue ξ.
4. Report `π̂ = mean(λ₀ + λ·h(y_i))` with interval `π̂ ± (β₀ + β_off + z·ŝ/√n)`.

- **Known-pose instance:** `ψ = Re<y − m_B(R), g(R)>/|g(R)|²` is unbiased under any zero-mean noise, with variance `(π(1−π) + E_ν[1/|g|²])/n`.
  - The price of view-law robustness is a harmonic, not arithmetic, mean of J.
  - On the three candidates `E|g|²` is 0.82, 1.50 and 0.28, giving standard errors of roughly 0.01–0.02 at 10⁴ particles with amplitude known *(est.)*.
- **Fit to the measured failures:**
  - integration error and likelihood tails affect efficiency only;
  - the viewing law is handled exactly;
  - amplitude is a latent coordinate;
  - a large required β₀ reports non-identification honestly.
- **Classical alternative:** your fixed-support joint mixture with a convex profile for π. It is the comparator, not the proposal:
  - its validity runs through likelihood values, grid support and boundary asymptotics;
  - its universal-inference version is certified only inside the finite-support model.

**E4. Theorem targets.**
- **T1.** Coverage uniform over μ, given a correct response. This is classical composition (Donoho; Armstrong–Kolesár).
- **T2.** The E2 identity and harmonic price with latent poses. This would be new, and I cannot promise it.
- **T3.** The identifiability curve, minimal variance against β₀, which is computable. The worst-case lower bound is nearly trivial: the adversary puts every particle on the least informative view.

**E5. Assumptions and exposure.**
- A closed two-state world, i.i.d. particles, and a forward model correct for first moments of bounded features.
- Off-catalogue bias can only be sampled. A rigorous Lipschitz bound needs a mesh near 0.04° *(est.)*.
- Weights scale like 4/J, so on 10076 the response must be accurate to a few 10⁻⁴ *(est.)*.
- Template error is not covered by any interval.

**E6. Novelty accounting.**

| Ingredient | Status |
|---|---|
| Per-particle marginal likelihood, ensemble weights | BioEM |
| Population intervals under a uniform view prior | cryo-BIFE (credible intervals) |
| Whole-stack population inverse problem | Evans et al.; replayed |
| Viewing-law misspecification; estimated view density | Xu et al., for structure with a bandlimited law |
| Mixture NPMLE, duality gaps, split likelihood ratio | Lindsay; Zhang et al.; Wasserman et al. |
| Bias-aware affine inference for linear functionals | Classical; empirical-Bayes and label-shift analogues exist *(not in packet)* |
| State-dependent views as unrestricted nuisance, with uniform coverage | Possibly new; the project's own reading does not establish absence |
| E2 as a diagnosis of pooled reweighting | Textbook with known poses; new only with latent poses |

- The statistics are a composition. ICML-level novelty would have to come from T2 or from a validated empirical finding, and neither is assured.
- The replay shows the released estimate moving from 0.80 to 0.65 under the condition labelled rotation-10. That shows pose-model sensitivity, not a viewing-law effect.

**E7. Computational bottleneck.**
- At 220 frequencies the template products are about 2×10¹¹ operations per stack, which is seconds. The response operator takes hours at most.
- The real limits are the off-catalogue check, in-plane shifts and higher resolution. Those confine a laptop study to domain-scale state differences.
- Splats help only as an exact off-grid template generator.

## 3. Protocol: three gates

Thresholds below are proposals; freeze your own before computing.

**G1. Materiality and identifiability (existing arrays and metadata; three days).**
- **Compute per stack:**
  - `J(R)` and amplitude-orthogonal `J⊥(R)` over the 65,536 saved rotations;
  - the odds factor `E_νA[J]/E_νB[J]` for a pre-listed set of recorded law pairs;
  - `E_ν[1/J⊥]`;
  - per-class view histograms on 10076 against the split-half floor of 0.025–0.03.
- **Pass, on at least two stacks:** predicted pooled bias of at least 0.02 at π = 0.75 for a measured pair, and `E_ν[1/J⊥] ≤ 9`, which is a standard error of 0.03 at 10⁴ particles.
- **My written prediction:** materiality fails for the compact regions with amplitude known. The moment ledger's view share of variance was only 15%, 3% and 0.5%.
- **Stop:**
  - if materiality fails, pooled likelihood is adequate and there is no methods contribution;
  - if identifiability fails, report occupancy at this band as unanswerable.

**G2. Frozen simulation under measured nuisances (laptop; two weeks).**
- **Design:** one simulator upgrade, frozen in advance, adding per-particle defocus, amplitude spread, shifts and coloured noise with whitening error. At most six scenarios, with 200 stacks of 8,192 each.
- **Comparators,** on identical images, templates, noise model and bank, tuned on separate seeds:
  - pooled reweighting;
  - a uniform-prior credible interval;
  - a shared estimated view density;
  - the joint fixed-support profile interval.
- **Pass:**
  - coverage of at least 185/200 in every scenario;
  - median width at most twice the known-pose oracle under Haar;
  - a pooled comparator below 0.80 coverage under state-dependent laws;
  - not dominated by the profile interval.
- **Stop:** end after one targeted revision fails. If the profile interval dominates, drop the methods claim.

**G3. Experimental contact (held-out particles).**
- **Controls on all three stacks:**
  - a phantom region in solvent, with upper limit at most 0.10;
  - a subtraction titration at 0, 0.25 and 0.5, judged on differences because baseline occupancy is unknown.
- **Truth:** one physically constructed mixture with source labels, covered at two of three mixing fractions, with the pooled comparator alongside.
- **Stop:** if controls fail on two stacks after one forward-model revision, make no experimental claim.

**Missing observations.**
- Noise-only micrograph windows for the in-band spectrum.
- A constructed mixture such as CAHRA challenge 1, including its frame-overlap details.
- Published class labels for 10076, if not on disk.
- No external compute is needed for G1 or G2.

**Unanswerable from current data.**
- The true viewing laws.
- Occupancy truth on 10028, 10049 and 10076. Three experimental stacks with population truth do not exist in the workspace.
- Template error.
- Solution-state populations, because of particle selection.
- A sufficient integration certificate on real images.

## 4. The catalogue revision

**Autopsy first, from saved arrays (free).** For each failing image:
- the share of total weight on Haar-labelled draws;
- the delta-method standard error against the between-bank difference;
- the planted-truth jump (N5);
- clipped eigenvalues (N7).

**If coverage or mis-weighting is confirmed, run it once.**
- Keep the local modes as the majority share. The 4.8° catalogue spacing cannot see a narrow secondary mode below about 5% of the dominant one unless weights are tempered *(est.)*.
- Fix one temperature, kernel width and share for all stacks before outcomes, and fit to the midpoint map only.
- The bank scores already exist from start selection.
- Write the prediction as lower-decile ESS of at least 900 (10049) and 450 (10076). 10028 should stay near 0.58 times the local share and still pass.
- Add the planted-truth report without changing the gate.

**Verdict.** Worth one day, in parallel with G1. It is classical repair and does not gate the scientific decision.
- **Pass:** competent marginal likelihoods for comparators.
- **Fail:** the branch ends, and comparators use fixed-bank scores with a stated limitation.

## 5. Priority order

1. Run the autopsy and G1 this week; both use existing arrays.
2. Add the D3 rows to the ledger, because amplitude-orthogonal information decides whether occupancy is viable at all.
3. Spend the single integration revision only if the autopsy supports it.
4. Build E3 only if G1 passes; otherwise write up the diagnosis and the negative result.
5. Request noise windows and a constructed mixture now, since G3 cannot start without them.
