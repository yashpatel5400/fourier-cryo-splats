I audited the enclosed source text without running anything. Findings below.

## Identities that check out

- **Derivative-tensor identity.** With `f = Re Σ c_q e_q`, write `Re A = (A + Ā)/2` and expand `∫‖D^m f‖²_F`. The cross pairing gives `<k_q^{⊗m}, k_l^{⊗m}>_F = (k_q·k_l)^m`, the phase factors give `(2πi)^m(−2πi)^m = (2π)^{2m}` for the difference term and `(2πi)^{2m} = (−1)^m(2π)^{2m}` for the sum term. The four blocks pair into conjugates, so the quarter becomes a half of a real part. The identity in the proposal and in `ball_sobolev_norms` is correct, and the sign `(−1)^m` on the sum-frequency term is right. Both kernels are real because both domains are origin-symmetric, which the identity needs.
- **Ball kernel.** The indicator transform is `4πr³(sin z − z cos z)/z³ = 4πr³ j_1(z)/z`. I expanded `sin z − z cos z` by hand through `z¹¹`; the five series coefficients in the code are correct and the dropped term is about `z¹⁰/5.2e8`, negligible below `z = 0.05`.
- **Cube kernel.** `Π_j sin(2πr v_j)/(π v_j)` equals `Π_j 2r·sinc(2r v_j)` with numpy's normalized sinc. Correct.
- **Path derivatives.** For `y(t) = R(t a u_r) x + t b` the derivatives are `a[u]_× R x + b`, `a²[u]_×² R x`, `a³[u]_×³ R x`. With a unit joint five-vector, `‖y'‖ ≤ aR‖u_r‖ + s_i‖u_t‖ ≤ hypot(aR, s_i)` by Cauchy–Schwarz, and `H = a²R`, `T = a³R` follow from `‖u_r‖ ≤ 1`. See E1 for the condition on the pose set.
- **Chain rule and change of variables.** `g''' = D³f[y',y',y'] + 3D²f[y'',y'] + Df[y''']`, pointwise Frobenius bounds, then Minkowski in L2 and a determinant-one rigid map onto a superset with nonnegative integrand. Sound. The Taylor factor `∫(1−t)²/2 dt = 1/6` is right.
- **Enclosures.** Ball: `‖Rx + tb‖ ≤ R + s_i`. Cube: `|(Rx + tb)_j| ≤ ½ + 2 sin(θ/2)‖x‖ + s_i` with `θ ≤ a`, monotone up to `π`, so the `min(angle, π)` cap is right for any angle.
- **Dual translation.** `k·b = q·τ` needs `E b = τ`, so `b = E⁺τ` and `‖b‖ ≤ shift/σ_min(E)`. The code's `lstsq`, `pinv`, spectral norm and the residual `2π‖q − kE⁺‖·shift` implement exactly this.
- **Phase-residual bound.** Leibniz on `(e^{itδ} − 1)M(t)` with `|e^{itδ} − 1| ≤ e` for `|t| ≤ 1` and `|M'| ≤ v`, `|M''| ≤ v² + h`, `|M'''| ≤ v³ + 3vh + j` gives exactly the recorded expression. Correct.
- **Pose convention.** `pose_derivative_fields`, `perturbed_geometry` and the test all use row `k @ exp(angle[u]_×)`, and the linearization matrix entries match `d/du_a (k·(e_a × x))`. The Hessian symmetric vectorization with `√2` is consumed correctly in the test.

## Issues

**E1. Pose-set mismatch between the two per-particle bounds. Severity: high, conditional.**
Location: `particle_ball_remainder`, default `joint_ball=True`; probe never passes it.
The old bound in `sharp_cube_cubic_coefficients` uses the sum speed `a‖k‖c + s‖q‖`, valid for a product pose set, and the docstring of `scaled_pose_radius` says the cubic remainder covers both sets. The ball bound uses `hypot(aR, s_i)`, which is only valid when `‖u_r‖² + ‖u_t‖² ≤ 1`. If the source audit was run with `pose_set='product'`, the per-particle minimum mixes a valid bound with an unproven one, and the recorded improvement is not a bound on the audited class. Nothing in the packet shows which set the source audit used.
Fix: pass `joint_ball = (audit['config']['pose_set'] == 'joint')`, raise if the key is absent, record the flag in the result JSON. Test: a product-set case with `‖u_r‖ = ‖u_t‖ = 1` where the direct remainder is compared against the sum-speed bound, and a regression test that the joint bound is strictly below the sum bound so the flag actually matters.

**E2. Radii hard-coded and additivity assumed. Severity: medium.**
Location: probe, `B, P = 2., 1.` and `bias = bias_upper − cubic_bias + selected`.
The 1e-12 check on the old cubic sum only verifies `B + P = 3`, not `B` alone, and `no_data = B·target_norm` depends on `B` individually. The subtraction is exact only if the source `bias_upper` is `min(joint, triangle)` with `(B+P)r` added to both branches, as in `joint_density_pose_bias`. The audit's `bound` dict is not shown, so this is assumed.
Fix: read `B`, `P` from the audit config; recompute `bias_upper` from the stored components with the old and new remainder and assert equality before use.

**E3. Missing targeted falsification tests. Severity: medium.**
Location: `tests/test_uq_ball_remainder.py`.
- Only eight random directions per case at one small angle. A bound that is loose by a large factor passes. Add a near-tight case: rotation zero, one frequency, real coefficient, where `g''' = −θ³ sin(·)` and the remainder norm has a closed form, and assert the ratio bound over actual lies in a known band.
- No test exercises `s_i = s/σ_min(E) > s`. Both nonlinear cases use an orthonormal embedding, so a missing `‖E⁺‖` factor would not be caught. Add a case with `E` scaled by about `0.3` and a large shift.
- No test above `angle = π` or near it, so the sine cap is untested.
- No test for `joint_ball=False`.
- The sum-frequency sign is only checked incidentally. Add a two-frequency case with `k_l = −k_q`, where `K(k_q + k_l) = K(0)` makes a wrong sign fail by a large margin.
- Add a coarse maximization over directions rather than random draws.

**E4. Test independence is partial. Severity: low.**
The nonlinear test subtracts `pose_derivative_fields(backend='direct')`, project code. The direct exponentials for the actual field are independent, so a wrong polynomial field would usually inflate the residual and fail, but this is indirect. A finite-difference check of the first and second fields along random directions would close the gap.

**E5. Floating-point guards are heuristic. Severity: low, disclosed.**
Location: `ball_sobolev_norms` pad `128·eps·N·mass`.
The pad ignores kernel evaluation error in `spherical_jn` and the relative error amplification in `dot**m` for `m = 3` at high band, where `(2π)^6‖k‖^6` is large. Recorded relative pads near 2e-11 mean cancellation is mild, so this is conservatism accounting, not a bug. A spot check of a few particles' `S_3²` in extended precision would substantiate the guard.

**E6. Forward-model sign convention unverified. Severity: low.**
`pose_cell_forward` applies `exp(−2πi τ)` while the adjoint fields use `+`. This is consistent only if `cell_forward` uses `exp(−2πi k·x)`. That module is not in the packet.

**E7. Provenance gaps. Severity: low.**
The ball-probe JSON lacks `domain` and neither result records `joint_ball` or the per-particle `path_speed_bound`. The proposal's first section still states `L = aR + s_i` and is only corrected later. The `isinstance(m, int)` check rejects numpy integers.

## Selection probability and event reuse

The per-particle minimum introduces no unaccounted selection probability. Both bounds are deterministic functions of the fixed weights, geometry, transfer and known noise scale. They do not depend on inference-set noise. The minimum of two valid deterministic upper bounds on the same quantity is a valid bound, and choosing to run the cube probe after seeing the ball result changes nothing for coverage. The noise event and the residual-norm event used by the source audit are untouched, so reusing them is legitimate. This holds only if E1 is resolved so both bounds cover the same pose class, and only under the source audit's assumption that weights were fit independently of the inference half.

## What the probes show

Both retained results are consistent with the code as read. The cube enclosure gives a smaller half-width than the no-data fallback, but the half-width still exceeds the target by more than an order of magnitude, and every reference scenario has zero correct-sign probability. Valid mathematics here changes a remainder constant; it does not change the density class, the assumed nuisance radii, or the absence of experimental calibration.

## Unverified

- The source audit's pose set and its `bound` decomposition.
- `cell_forward` sign convention.
- Accuracy of `spherical_jn` and pairwise summation at the recorded magnitudes.
- Weight independence from the inference half.
- The test log itself, which I did not run.
