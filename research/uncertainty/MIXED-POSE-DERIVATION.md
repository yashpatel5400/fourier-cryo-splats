# Conditional mixed pose-error bound

1 October 2026 UTC, response to full review 2 (R9). This is a new development lemma assembled from the existing continuous derivative fields and the classical Hoeffding bound. It is not a claim of novel concentration theory or calibrated experimental poses. It has numerical tests, but has not yet received external review.

Let the nominal design, weights w, density rho=rho0+h, ||h||≤B, and pilot be fixed. The pose error for particle i is U_i=mu_i+xi_i. Assume the xi_i are independent, centered, independent of Gaussian measurement noise, and ||xi_i||≤r_i. Allow arbitrary common/dependent components mu_i with ||mu_i||≤eta_i. All these statements must hold conditional on whatever was used to choose the design and weights. They do not follow from local alignment on the same images.

Write the error of an affine feature estimate as its nominal density residual, plus first-order pose terms, Taylor remainder, and measurement noise. Let D_i:R^5→L2 be the derivative field built from particle i's estimator weights, and v_i=D_i* rho. A valid scalar bound is

c_i = ||D_i* rho0|| + B ||D_i||op,

or c_i=(||rho0||+B)||D_i||op without pilot pairings. The common-mode contribution is bounded by sum_i eta_i c_i. Hoeffding's lemma gives E exp(t v_iᵀxi_i) ≤ exp(t² r_i² c_i²/2). Multiplying these mgfs and the independent measurement Gaussian mgf gives variance proxy s²+sum_i r_i² c_i², where s²=||w||².

For nonlinear rotations/translations, Taylor's formula gives a remainder bound

R2 = (||rho0||+B)/2 sum_i (r_i+eta_i)² sum_q |c_iq| [(2pi(angle |k_iq| (13/720)^(1/4) + shift |q_iq|))² + 2pi angle² |k_iq|/sqrt(12)],

where c_iq=(w_real+i w_imag)CTF/noise. The directional fourth moment of a unit-cube coordinate is at most 13/720, and the second is 1/12. Rotation preserves |k| along the Taylor path, so this bound is uniform over each full pose segment; it is not only a derivative evaluated at the nominal pose.

Consequently, with probability at least 1−alpha,

|estimate−true feature| ≤ B||ell−A*w|| + sum_i eta_i c_i + R2 + sqrt(2 log(2/alpha) [s²+sum_i r_i²c_i²]).

This uses a sub-Gaussian tail, not the standard normal quantile suggested informally in the review. No Gaussian assumption is made on xi_i. First derivative Gram entries are analytic moments of the continuous cube using derivatives of sinc. The implementation reports ordinary floating-point diagnostic padding, not a validated rounding enclosure.

Under repeated identical designs and weights scaled by 1/n, the stochastic first-order proxy and measurement variance scale as 1/n. The common-mode term and quadratic remainder need not shrink. Thus this change alone does not remove systematic alignment error or nonlinear bias, and cannot promise asymptotic vanishing widths.

## Estimated-pose limitation to test

For weights/design computed from an alignment image, the pose errors generally are neither centered nor independent of those weights conditional on that image. Splitting alignment and inference noise fixes the measurement-noise dependence, but by itself does NOT prove the centered-pose premise. The forthcoming simulation will explicitly report empirical coverage of this plug-in use, including failures; it cannot promote this conditional lemma to an end-to-end theorem merely because coverage looks favorable on a few generators.
