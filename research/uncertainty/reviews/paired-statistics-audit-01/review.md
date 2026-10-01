# Focused audit: paired-exposure power and matrix-covariance candidates

**Bottom line.** The Gaussian algebra, cone lemmas, translation claim, `svec`/`smat` and LP dual signs all check out on paper. The diagonal negative result on 10076 stands. The matrix screen as declared cannot deliver a decisive answer in either direction: it has no upper bound, an unfavourable basis, no held-out-view check and a fixed translation.

## Identities that check out

- **Scalar/complex bilinear moment.** Conditioning on Z gives `E[exp(tXZ)|Z] = exp(tmZ + t²vZ²/2)`. Integrating over Z gives exponent `[t²vm²/2 + tm² + t²vm²/2]/(1−t²v²) = tm²/(1−tv)` and prefactor `(1−t²v²)^(-1/2)`. Two real coordinates give `log(1−t²v²)` and `u|μ|²`, and the translation phase cancels.
- **Covariance monotonicity, signed T, noncommuting S.** `det(I−2HS) = det(I − T S_Z T S_X)`, and with `A = T S_Z T ⪰ 0` this is `det(I − A^(1/2) S_X A^(1/2))`, decreasing in `S_X`. `K = (I−2HS)^(-1)H = H(I−2SH)^(-1)` is symmetric and `dK = 2K(dS)K ⪰ 0`.
- **A shorter proof of the same fact.** `E[exp(X'TZ)|Z] = exp(m_X'TZ + ½Z'T S_X T Z)` is pointwise monotone in `S_X`; then swap roles. This avoids the path argument, covers unequal means, and should replace the derivative argument in the text.
- **Moment domain.** `‖S_X^(1/2) T S_Z^(1/2)‖² ≤ λ_max(D^(1/2) T D T D^(1/2)) = max(t_j v_j)² < 1`. This uses that diagonal T and D commute; in the matrix case `D = I`. The whole path stays in the domain.
- **Mismatch coefficient.** The exponent `[t Re(μ_x μ̄_z) + (t²v/2)(|μ_x|²+|μ_z|²)]/(1−t²v²)` matches the code. It reduces to `t/(1−tv)` at r=1, φ=0. It is convex in r, and the cosine choice for t<0 is right. Unconstrained phase makes every nonzero coefficient positive, as claimed.
- **Matrix form.** An orthogonal change of basis preserves `S = I`, so the scalar product formula gives `det(I−T²)^(-1/2) exp(m'Um)` with `U = T(I−T)^(-1)`. `matrix_gaussian_terms` implements this.
- **Cone lemmas.** `U−T = T²(I−T)^(-1) ⪰ 0`, so `tr(TS) ≤ tr(US) = Σλ_l m_l'Um_l ≤ 0`. The diagonal analogue `t/(1−tv) − t = t²v/(1−tv) ≥ 0` holds. These are statements about expected log growth only.
- **Small-c separator.** `T = U(I+U)^(-1)` inverts `U = T(I−T)^(-1)`. The domain `λ_min(U) > −1/2` and the derivative `tr(HS)` at c=0 are correct.
- **Translation claim.** Commutation with all block rotations kills the complex-linear part of block (j,k) unless `q_j = q_k`, and the antilinear part unless `q_j = −q_k`. A diagonal block is `aI + bJ`, and symmetry forces b=0. Bounded shift sets suffice, by analyticity. The dropped `J` term is `Im(X Z̄)`, which has zero mean under common means.
- **`svec`/`smat`.** The √2 scaling preserves the Frobenius inner product, and the inverse is consistent.
- **LP dual sign.** SciPy/HiGHS marginals for `A_ub x ≤ b` in a minimisation are ≤ 0. With `y = −marginals`, dual feasibility gives `a_l·(y₂−y₁) ≤ 0` per view and objective `b·(y₂−y₁) = e*`. The code's `(marg[:n] − marg[n:])/scale` equals `(y₂−y₁)/scale`, so `⟨H, m_l m_l'⟩ ≤ 0` and `tr(HS) = e* ≥ 0`.
- **Repairs.** Subtracting `c·I` on trace-normalised views (matrix) or a constant on row-normalised powers (diagonal) shifts every constraint by exactly c. The fallback catches `u ≤ −1/2`.
- **Diagonal program and dual.** `x/(1−x) = inv_pos(1−x) − 1` is convex and the objective is concave. The dual is separable, and any nonnegative multiplier gives a weak-duality upper bound.

## Findings

**P1 — High — `covariance_cone_diagnostic` + `direction_growth`; MATRIX-PROPOSAL "Declared first screen".** The matrix screen reports only a feasible lower value along one direction. That direction is the dual of a coordinate-scaled ℓ∞ LP, which optimises an ℓ1-normalised margin, not growth. Unlike the diagonal gate there is no upper bound, so small growth cannot support a stop and near-membership is unquantified.

*Fix:*
- Replace the LP with a Frobenius projection: `A* = argmin ‖S − Σλ m m'‖_F`, λ ≥ 0 (NNLS in `svec` coordinates).
- By Moreau decomposition, `Δ = S − A*` is itself a separator: `⟨Δ, m_l m_l'⟩ ≤ 0` and `⟨Δ, S⟩ = ‖Δ‖²`.
- It also gives a closed-form upper bound. For any feasible T, `tr(TA*) ≤ tr(UA*) ≤ 0`, so growth `≤ sup_T tr(TΔ/v) + ½ log det(I−T²) = Σ_i g(δ_i/v)`. Here `g(δ) = δt* + ½log(1−t*²)` with `t* = (√(1+4δ²)−1)/(2δ)`, and `g ≈ δ²/2`. Von Neumann's trace inequality makes the supremum commute with Δ.
- The exact problem is also convex (`m'(I−T)^(-1)m` is matrix-fractional and the log-det terms are concave), so ranks 8–16 can be solved directly.

**P2 — High — `probe_uq_paired_covariance.py` basis choice; MATRIX-PROPOSAL "top 32 eigenvectors".** The basis maximises shared signal energy, which is where the true and region-removed maps agree most. A local removal lives in low-eigenvalue directions, so `retained_true_energy` near 1 is compatible with zero retained discrimination.

There is also a structural mechanism. `S_alt` is the null orbit's uniform average (an interior point of the null cone, if the cone is full-dimensional in Sym_r) plus a small perturbation. At r = 8/16/32 the symmetric dimension is 36/136/528 against 4,160 views, so full dimensionality is likely and small alternatives are then automatically inside. This is the same mechanism that plausibly produced the 10076 diagonal witness.

*Fix / falsification:*
- Report the numerical rank of the feature matrix against r(r+1)/2.
- Add a discrimination-oriented oracle basis, for example leading eigenvectors of `S_true − S_removed` or a generalised eigenproblem.
- Report `‖P_B(S_true − S_removed)P_B‖_F / ‖S_true − S_removed‖_F`.
- Do not read a negative at these ranks with this basis as a negative for the matrix family.

**P3 — High — repair steps in `uq_power_cone.py` and `uq_paired_covariance.py`; both propositions.** Amplitude is unbounded, so the null expectation at an off-catalog view is `exp(a² m'Um)`. Any strictly positive off-grid violation makes the null expectation unbounded, not slightly above 1. The `1e-12` margin gives no slack between views. The separator has up to 528 free parameters against 4,160 constraints, and the diagonal gate had d parameters against 64 views. Positive finite-view growth is therefore partly catalog overfitting; the diagonal result already shrank this way.

*Falsification test:* evaluate each repaired U on at least 10⁴ fresh Haar rotations of the null candidate. Report the maximum `m'Um/‖m‖²` and the fraction positive. Recompute growth after repairing against the union, and iterate until growth stabilises or vanishes. This is cheap and should precede any SO(3) enclosure work.

**P4 — High — MATRIX-PROPOSAL "Translation tradeoff" / "initial screen fixes translations exactly".** The invariance claim is correct, but the fixed-shift screen gives the alternative perfect centering and the null no shift freedom. With a shift distribution, entry (j,k) of `S_alt` is damped by the characteristic function at `q_j − q_k`: `exp(−2π²σ²|Δq|²)` for Gaussian σ. Only pairs with `|Δq| ≲ 1/(2πσ)` keep phase information, and broad shifts collapse `S_alt` to the already-defeated diagonal. The null cone also grows to a 5-parameter orbit.

*Gate:* rerun the P1 bound pair with a declared alternative shift law (σ at 0.5, 1, 2 px) and a null catalog covering rotations × bounded shifts. If growth does not survive realistic σ, stop.

**P5 — Medium — FINITE-VIEW-RESULTS closing paragraph; ENLARGED-CONE-RESULTS "reason to stop".** The packet separates expectation from power correctly, but three sharpenings are needed.
- Cone membership gives nonpositive conditional expected log for every admissible predictable T. Log-wealth is then a supermartingale under the alternative, and a fixed `T ≠ 0` i.i.d. product tends to 0 almost surely. This is stronger than "no positive expectation". It still does not bound finite-n rejection probability: `E_alt[E]` can exceed 1 by Jensen, so Markov does not apply. Nor does it rule out averaged or mixture e-values.
- "Below log(20) at 128 particles" is not a meaningful criterion, because growth is linear in n and real stacks have 10⁴–10⁵ particles. The decisive evidence is membership, and it exists only for 10076.
- 10028 is trending toward membership (.0425 → .0088) but is unproven. 10049 (.288) is not defeated by this packet. State the stop decision as "defeated on 10076, plausibly on 10028, open on 10049".

The 10076 conclusion does not need an algebraic identity. The dual value (7.55e-28) is the bound, and it is quadratic in the residual.

**P6 — Medium — `cone_witness` (zero objective).** The reported mass 1.0926 is whatever vertex HiGHS returned, and the implied 4.5% amplitude excess is not canonical. *Fix:* solve min-mass and max-mass LPs subject to the equality. If min mass ≤ 1, the failure survives a null with amplitude bounded by 1; otherwise the result quantifies exactly how much amplitude freedom the failure needs.

**P7 — Medium — `covariance_cone_diagnostic` dual extraction.** The sign is correct, but nothing asserts it at run time. An inaccurate or mis-signed dual produces a large `correction`, an almost negative-definite direction, and zero growth, which looks the same as true non-separation. When the target is in the cone the dual is degenerate, and `raw_direction/max(norm, 1e-30)` can normalise roundoff to a unit-norm arbitrary direction. That is harmless by the lemma but meaningless.

*Fix:* assert `|dual_objective − lp_objective| ≤ tol` and `raw_maximum_normalized_violation ≤ tol`. Flag cells where `correction·tr(S)` is comparable to `e*`. Return a zero direction when `lp_objective < tol`.

**P8 — Medium — `probe_uq_paired_covariance.py` versus the stated model.** The script is a moment computation; no X or Z is drawn. Expected log growth does not depend on the true noise covariance, so the Loewner-bound machinery is exercised only by the 5-dimensional unit test. Specific gaps:
- "First saved CTF/noise profile" applies only `transfer_squared[0]`. No noise spectrum or whitening appears, so identity covariance in these coordinates is assumed, not implemented.
- If the saved frequency set contains Hermitian pairs (q and −q), the realified covariance has eigenvalue 2v and violates `S ≤ D`. This would invalidate the diagonal proposition too. The packet does not establish nonredundancy.
- Self-conjugate frequencies have one real coordinate. The complex normaliser is then conservative, but this is undocumented.

*Test:* assert nonredundancy on the frequency arrays. Run a Monte Carlo null check of `E[E] ≤ 1` at catalog and held-out views with noncommuting covariances at the bound.

**P9 — Low/Medium — `scalar_dual_maximum`.** The value is evaluated at an approximate root, so it is at most the true supremum. The labelled "upper" is slightly anti-conservative. This is negligible in the interior but uncontrolled as x → 1 (large s, zero-signal rows). *Fix:* add `|f'(x̂)|·bracket_width`, or evaluate a tangent bound from both bracket ends.

**P10 — Low — tests and replay.**
- `sign > 0` in `dense_mgf` is not a domain check; an even number of negative eigenvalues would pass.
- No matrix test covers unequal means, eigenvalues near ±1, or singular covariances.
- The translation test covers one frequency's diagonal block, not the off-diagonal claim.
- `verify_uq_power_cone_witnesses.py` reuses `power_witness_log_upper`, so the replay is independent of the LP but not of the dual evaluator.
- `direction_growth` is not concave in c. Its values are valid lower bounds regardless, but "locally refined" should not be read as the directional maximum.

## Most consequential missing gates, in order

1. **Held-out-view violation (P3).** Minutes of compute. It decides whether any finite-view positive is real.
2. **Two-sided bound with a discrimination-aware basis (P1, P2).** This turns the matrix screen into a stop/go test.
3. **Translation-distribution rerun (P4).** This decides whether the matrix extension has any content beyond the diagonal.
4. **Power at realistic n.** Compute `Var(log E)` analytically for the selected T, apply a normal approximation at n = 10³–10⁵, and run one Monte Carlo product check. Do this in place of further expectation tables at n = 128.
5. **Null calibration with the true map.** Use adversarial covariances and both catalog and off-catalog views.
6. **Paired empty-ice cross-power.** Solvent, neighbouring particles and fixed-pattern errors are common to both exposures. They add positive cross-power, which would falsely reject a true map through positive weights. Zero-mean cross-power on particle-free windows is a cheap test of the independence premise.
7. **An oracle comparator.** Run a known-pose likelihood-ratio or half-map test on the same alternative at the same n. If the oracle bilinear family cannot approach it, stop without a continuous certificate.

## What was and was not verified

**Verified by hand from the packet text:** everything under "Identities that check out", and the code-to-formula correspondence in the five source files as printed.

**Not verified:**
- No code was executed, and no arrays, logs beyond the printed text, or original publications were inspected.
- I did not check the frequency-set nonredundancy, the noise whitening, the number of frequencies or the catalog construction.
- I did not check the cvxpy or HiGHS numerical behaviour, the 288-case and 36-case numbers, or the six tests' actual pass status.
- The matrix screen has no results in this packet, and I assert none.
- Nothing here establishes continuous-orientation validity, rejection power, experimental calibration or novelty.
- The negative conclusion for the diagonal statistic on 10076 is preserved.
