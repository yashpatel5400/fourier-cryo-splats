# Candidate: paired-exposure Fourier-power compatibility tests

1 October 2026 UTC. Method development after the frozen review-3 snapshot; no dataset-level power or coverage experiment has been run for this candidate. This changes the target to compatibility of a candidate map's projection powers, not a confidence interval for a local density. It does not repair the failed v0.7 pose bounds. Novelty is unestablished; Gaussian moment-generating functions, e-values, moment-based cryo-EM comparison, alignment-error models and independent-particle validation are prior art.

## Motivation and scope

The local-pose study exposes large nonlinear remainders and uninformative joint pose balls. Instead of estimating each pose and bounding its error, retain every possible orientation in the null. Fourier power removes translation phase. Allow unknown nonnegative particle amplitude. Separate repeated exposures offer a bilinear statistic whose Gaussian noise normalization can use a covariance upper bound; an estimated pose radius and a density-deviation norm radius are absent.

The price is a weaker estimand. Failure to reject does not establish the density, local feature accuracy, a unique structure or phase recovery. Translations, global rotations, handedness and power-equivalent alternatives must be negative controls. An exact candidate-map null will generally be false for an approximate reconstruction; a rejection is incompatibility under the observation model, not a biological truth label.

## Conditional model

Condition on an independent training set that fixes the candidate Fourier map F, frequencies, weights, processing and metadata. For held-out particle i, take two realified complex Fourier observations X and Z with the same mean mu and independent Gaussian errors. Within either observation, arbitrary correlations are allowed provided both covariance matrices are bounded in Loewner order by D = diag(v_1,...,v_d). For the convenient complex formula below, the two real coordinates of frequency j share upper bound v_j. These are upper bounds on the full covariance, not just its diagonal entries.

For candidate orientation R, shift s and arbitrary amplitude a >= 0,
mu_j = a C_j exp(-2 pi i q_j.s) F(R q_j).
Repeated-exposure means must agree after any independently specified dose/gain corrections. Independence of their errors and independence of the training/processing graph remain substantive requirements. The existing raw-movie diagnostic does not establish them.

Choose real betting weights t_j with |t_j| v_j < 1. Define u_j = t_j/(1-t_j v_j). A sufficient candidate-map constraint is

    sup_R sum_j u_j |C_j|^2 |F(R q_j)|^2 <= 0.

It is a constraint on a map's continuous orientation orbit. Satisfying it only on a sampled rotation grid is not a continuous certificate. When it holds, define

    log E_i = sum_j [t_j Re(X_j conj(Z_j)) + log(1-t_j^2 v_j^2)].

### Proposition (candidate validity under the stated Gaussian model)

For every orientation, translation and nonnegative amplitude in the null, E[E_i | training, metadata] <= 1. Products over independent particles, or a predictable sequence satisfying the same conditional assumptions, give an e-value. Rejecting when the product is at least 1/alpha controls Type-I error by Markov's inequality. Learning the weights on the inference pixels is not allowed by this proposition.

### Derivation

For one real coordinate with equal variance upper bound v and common mean m, the exact bilinear Gaussian moment at the bounding covariance is

    E exp(t X Z) = (1-t^2 v^2)^(-1/2) exp[t m^2/(1-tv)].

For general independent within-image covariances S_X,S_Z <= D, write the quadratic form on (X,Z) with off-diagonal blocks T/2. Its central determinant factor is det(I-S_X T S_Z T)^(-1/2). This increases when either covariance increases in positive-semidefinite order: use the determinant identity to express each comparison through the positive-semidefinite matrix T S_Z T or T S_X T. The common upper covariance D consequently supplies an upper bound even for signed T.

The noncentral quadratic coefficient is K(S)=(I-2 H S)^(-1)H, which is symmetric on the valid moment domain. Along a positive-semidefinite covariance increment, its derivative is 2 K(S) (dS) K(S), also positive semidefinite. Thus the same bounding covariance upper-bounds the noncentral exponent. The domain |t_j|v_j<1 guarantees finiteness for all smaller block covariances. Combining the two real coordinates at each frequency gives the displayed log normalizer. Its remaining noncentral exponent is a^2 times the candidate-map constraint; translation phases vanish. This proves the proposition.

The elementary Gaussian identities and e-value conversion are not asserted as novel. The unresolved contribution would be an efficient, useful continuous-orientation constraint and an experimentally defensible paired-exposure pipeline.

## Relative frame motion and dose are still nuisances

A common translation cancels from the cross power; different frame translations do not. This is an essential qualification. For a second mean r exp(i phi) times the first mean, with r in [r_low,r_high] and |phi| <= epsilon, replace u at each frequency by

    max over r in {r_low,r_high} of
    [t r c + (t^2 v/2)(1+r^2)] / (1-t^2 v^2),

where c=1 for t>=0 and c=cos(min(epsilon,pi)) for t<0. The expression is convex in r, so endpoints suffice. With independently established bounds on relative gain and phase, substituting this upper coefficient in the orbit constraint preserves the proposition. It can be conservative because separate frequency maxima need not be jointly attainable.

An unconstrained relative phase makes every nonzero coefficient positive: the cancellation needed for a nontrivial amplitude-free power-cone witness disappears. The method therefore avoids particle viewing-angle estimates, but still needs an evidence-based relation between paired exposure means. This limitation must be included in any physical-data protocol.

## Primitive numerical checks, not study results

Six dedicated tests pass. They compare the scalar/complex formulas against an independent dense block-Gaussian determinant/inverse calculation, including 100 predetermined random cases with noncommuting within-frame covariances below the upper bound, signed weights, gain/phase mismatch, exact endpoint equality, zero-noise limits and the unconstrained-phase failure. They do not certify an orientation supremum, demonstrate statistical power or validate experimental assumptions. The historical 279-pass suite was not rerun and its reported count is unchanged.

## Falsification before large experiments

1. Independently verify the Gaussian algebra for noncommuting, correlated covariance matrices below D, signed weights, nonzero means and boundary cases. Compare with direct block-Gaussian determinants/inverses; retain failures. This is a numerical implementation check, not coverage validation.
2. An oracle finite-view separation probe can upper-bound potential discrimination: a sampled-view null is weaker than the continuous null. If even that favorable relaxation cannot separate declared structural alternatives, stop this candidate. A favorable relaxation does not certify a continuous test.
3. Only if justified, build a continuous SO(3) upper bound. Interval bounds on each Gaussian-splat Fourier power can provide sufficient constraints. Their looseness, runtime and numerical enclosure must be measured; no sampled maximum may be relabeled a supremum.
4. A new protocol must then freeze independent fitting/calibration/test seeds, null transformations, amplitude controls, corruption alternatives, all three geometries and comparison methods before test outcomes. The previous 600 trials cannot be reused as a fresh confirmation.
5. Experimental tests require independent exposure processing, defensible covariance bounds, and treatment of CTF/dose/window effects. No claim of solving those inputs follows from the proposition.

Closest cryo-EM moment-comparison and likelihood-validation sources are already discussed in `../MOMENTS-AND-IDENTIFIABILITY-READING.md` and `../LIKELIHOOD-VALIDATION-READING.md`. The additional alignment-moment reading is in `../ALIGNMENT-MOMENT-FOLLOWUP.md`. The proposal and all subsequent results must distinguish a map-compatibility test from density coverage.
