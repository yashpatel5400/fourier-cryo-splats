# Follow-up audit: pilot pairing refinement, B9 to B11 fixes, and accounting

No tools were available; every claim below was re-derived by hand from the packet. I found no invalid upper or lower bound. Three low-severity issues are new.

## Pilot pairing derivation

The scalar bias is `<u,h-Fv-e> - <rho0,Fv+e>`. The only changed term is `<rho0,Fv>`. Both `|<rho0,Fv>| <= P||Fv|| <= P f` and `|<rho0,Fv>| <= ||F*rho0|| ||v|| <= p L` are Cauchy–Schwarz, so `min(P f, L p)` is valid, and the density interaction, cross term and remainder are untouched. The density class is unchanged because `u` never enters the modified term. Since `p` is deterministic, taking the minimum with the spectral-event bound `f` and then with the triangle form stays valid on the same event, which is what `joint_density_pose_bias` does (`pilot=min(P*f, L*pilot_cross_upper)` applied to both `joint` and `triangle`).

The scaling is consistent: `pair_pose_fourier_moments` divides by `op.denominators`, so `p` is the adjoint norm in the same lifted coordinates in which `||v|| <= L = sqrt(sum d)` and in which `c = ||F*h||` is computed. The contraction (constant/linear phase coefficients for first-derivative columns, `-product + i curvature` with `.5*scale` for the fifteen half-Hessian columns) is byte-for-byte the `rmatvec` contraction with continuous moments substituted for quadrature sums, and the existing Gaussian test compares exactly that substitution against `rmatvec`.

**Archived numbers check.** For the adaptive-average pilot audit: `L = sqrt(19.0324) = 4.3626`, `L p = 4.3626 × 0.018144 = 0.07916`, matching `new_pilot_polynomial_upper`. Bias `8.61353 - 0.83491 + 0.07916 = 7.85778` matches. Half-width `7.85778 + 1.6449 × 0.24104 = 8.2543`, relative `8.2543/16.180 = 0.5102`, change `8.2543/10.447 = 0.7901`. The stress ratio `0.6667 × 7.8578 = 5.239` equals the earlier `0.5818 × 9.004` for the same saved poses. All consistent.

## Constant-cell moments, NUFFT signs, half-cell phase

With `rho0 = m^{3/2} sum_cells coef 1_cell`, and `x = c + z` on a cell of volume `m^{-3}`:

`∫ rho0 x^β e^{2πik·x} = m^{-3/2} sum_cells coef e^{2πik·c} ∏_a sum_g C(β_a,g) c_a^g E[z^{β_a-g} e^{2πik_a z}]`.

The code implements exactly this: the binomial factor is only nontrivial for `β_a=2, g=1` (the explicit 2), the `m^{-3/2}` appears as `/box**1.5`, and cell centers `(i+1/2)/m - 1/2 = (i - m/2)/m + 1/(2m)` give `e^{2πik·c} = e^{i(2πk/m)·mode} e^{iπ k/m}`. Type-2 FINUFFT with `isign=1`, nonuniform points `2πk/m`, and default mode ordering from `-m/2` reproduces the first factor, and `exp(1jπ k.sum/box)` is the half-cell phase with the correct sign. The (z,y,x) array order is paired with frequency axes `[2,1,0]`, which is self-consistent under the FINUFFT Python convention that array axis 0 corresponds to the first coordinate; a wrong pairing would be caught by the composite-cell test with random asymmetric coefficients.

**Scalar formulas.** With `x = πk/m` and NumPy `sinc(u) = sin(πu)/(πu)`: `E[e^{2πikz}] = sinc(k/m)`, `E[z e] = sinc'/(2πi m)`, `E[z² e] = -sinc''/(4π² m²)`. `sinc'(u) = π(x cos x - sin x)/x²` and `sinc''(u) = π²(-sin x/x - 2cos x/x² + 2 sin x/x³)` match the code. I expanded both series: `-1/3, 1/30, -1/840, 1/45360, -1/3991680` and `-1/3, 1/10, -1/168, 1/6480, -1/443520` are correct, and the `|x| < 0.03` switch with terms through `x^{10}` leaves truncation far below `1e-15`. The scalar test places points at 0, `1e-8`, `-1e-5` and on both sides of the switch.

**Independence of tests for normalization.** The composite test integrates `coef · m^{3/2} · x^β e^{2πik·x}` with a 12-node rule inside every cell, using its own centers, its own `m^{3/2}` and its own weights (`1/m³` per cell). A wrong `box**1.5`, a wrong half-cell sign, a missing binomial 2, or a swapped axis all fail it at `rtol 2e-9`. The zeroth moment is additionally checked against `cell_forward`'s independent conjugate convention. The finite-difference test checks the pairing against `pose_cell_forward` for random pilots: the first difference must equal `<rho0, D_a'w>` and the second difference must equal `uᵀHu`, which requires the `.5*scale` and `svec` conventions to be right (I verified `2·pair@svec(uu') = uᵀHu`). These three checks share no normalization code with `cell_fourier_moments`. The Hilbert-space test only checks the inequality, not normalization, which is appropriate.

## B9 to B11 verification

**B9 (closed).** `recover_reference_centers` now asserts both broadcast radius arrays are constant and passes `float(...flat[0])` to `perturbed_geometry`, whose `angle*pose[:,:3]` and `shift*einsum(...)` then broadcast correctly. Nonuniform radii raise rather than silently pick one value. The regression fixture deletes `pilot_target` and `raw_expected_center` and forces `uses_no_data=True` on the archived 10049 two-degree diagnostic, so it takes the full reconstruction path with the same seed rule (`609711+10049`) as the diagnostics script; the recovered means equal the archived raw means to 0.0, as expected for a deterministic replay. The fixture passes `np.full(128, ...)` arrays, which is the same shape the audit's call site supplies. This closes the functional defect.

**B10 (closed).** `previous.get('pilot_target')` yields `None`, and `reference_interval_summary` raises on fallback with `None` while permitting it on the affine branch; the test covers both.

**B11 (mostly closed; see B13).** `audit_population` scans every JSON under the development tree for the audit stage string and lists the source-fit grid with unfinished fits captured by hash. The main table's ten-row scope is stated.

## Design-order scale selection (re-confirmed)

For any positive block scales `d`, `sum H_j v_j = F̃ u` with `F̃ = [H_j/sqrt(d_j)]` and `||u||² = sum d_j||v_j||² <= sum d_j`, so `||field|| <= sqrt(sum d) ||F̃||_op`. In `audit_uq_high_band_pose.py` the operator is built at `args.order`, `establish_quadrature_design_scaling` writes only `self.denominators`, the probes act on the final-order `spatial_gram`, the pad uses `polynomial_kernel_error(k, args.order)` with masses divided by the same denominators, and `L = scaled_pose_radius(group_scales)` uses the same scales. Sampled low-order norms never enter any bound. The guarantee is preserved; poor scales only loosen it. The pilot pairing in that script uses the same `op`, so `p` is in the same lifted coordinates.

## New issues

**B13. Packet summary is stale relative to the summary script.** Severity: low (reporting). The included `summarize_uq_joint_bias.py` emits a `pilot_pairing` key in `all_attempts`, but the included `summary.json` has no such key and lists 22 attempts, so it predates the pilot audit. The 0.5102 pilot-refined width appears in no inventory. Fix: regenerate the summary and state that pilot-pairing widths are a separate variant of the same conditional bound.

**B14. Pilot pairing norm carries no numerical pad.** Severity: low (disclosed numerics, not validity in real arithmetic). `p` is an exact continuous quantity, so no quadrature pad is owed, but its NUFFT evaluation (`eps=1e-12`) has no archive-scale independent check comparable to the stress test's direct-sum `independent_field_relative_error`. At `p = 0.018` and `L = 4.36`, a relative NUFFT error of `1e-10` changes the bias by `~1e-11`, so nothing archived is at risk. Targeted check: recompute a handful of lifted columns at archive scale by direct exponential sums over the `24³` centers and record the relative error in the audit JSON.

**B15. High-band audit assumes a pose-free source bias.** Severity: low. `fit['bias']/B` is used as `||h||`, which is correct only if the source is a fixed-pose fit whose bias is `B × padded residual`. The script checks completeness but not that the fit lacks pose terms. Fix: assert `'pose_polynomial_bias' not in fit` or read `density_bias` explicitly.

## Identities confirmed correct

- `min(Pf, Lp)` bound and its triangle counterpart, on the same spectral event.
- Constant-cell moment factorization, binomial expansion, `m^{-3/2}` normalization, mode ordering and half-cell phase sign.
- Scalar `sinc'`, `sinc''` closed forms and all ten series coefficients.
- Cross-error pad `E T ||m||` (unchanged), joint and product radii (unchanged).
- `2·pair@svec(uu') = uᵀHu` so the finite-difference second-derivative test is correctly scaled.
- B9 assertion and scalar conversion; B10 `None` semantics.

## What I could and could not verify

Verified from the packet: the pilot-pairing inequality and its implementation, the cell-moment derivation and code, near-zero series, the test constructions and their independence from the moment code's normalization, the archived pilot-audit arithmetic, the B9 fix and its regression fixture logic, the B10 change, the population scan logic, and the design-order argument.

Not verifiable here: any execution (the 13- and 93-test logs, the 0.0 regression difference), FINUFFT's axis convention and accuracy at `24³` modes, the `model(...)` helper and whether the pilot checkpoint is the one the source fits used, `pose_cell_forward`'s rotation sign agreeing with the polynomial operator's convention (the finite-difference test would fail otherwise, but I cannot run it), and file hashes. The refined widths remain conditional bounds for fixed weights on the stated joint ball, with an unpadded but tiny numerical pilot pairing; they are not feature detection, calibration, or experimental coverage.
