# Audit-motivated classical risk comparators and slack diagnosis

This post-outcome saved-array analysis addresses focused Fable audit V1--V5.
It is not new independent calibration or experimental validation. All four
ranged-amplitude scores, both candidates, all five kappa values and three
sample sizes are retained. Existing arrays and their hashes are fixed before
this script runs. No minimum across separately valid methods is used as a test.

## Assumption correction (V1)

For the grouped and variance procedures, physical amplitude A and measurement
noise epsilon must be conditionally independent given the view R. Amplitude
may vary with view and between particles, but entire tested particles must be
independent. Image-dependent selection or an image-fitted scale is not covered
merely by calling it an amplitude. The individual-noise envelope protects the
larger model in which amplitude can depend on that particle's noise, assuming
the stipulated noise marginal and independent particles still hold. Comparisons
between individual and grouped envelopes therefore also change this part of
the null model. The existing simulator does meet the stronger independence
assumption. This addendum corrects an omitted qualifier; it does not rewrite the
frozen experiment or claim the qualifier was experimentally established.

Write eta_L(R)=E[X|R], explicitly including finite L=32 and cell-max bias.
For any w=dP/dQ with E_Q w=1 and E_Q(w-1)^2 <= r, the variance-only branch
gives E_P eta_L <= mu+sqrt(r V). It does not require a pointwise density cap.
Under this weaker chi-squared assumption alone, the kappa*mu branch must be
discarded. A cap kappa implies r<=kappa-1, as in the archived study. No
experimental chi-squared radius is estimated or claimed easier to calibrate.

## Comparators (V3)

The exact bounded-density population optimum is the upper-tail mean
CVaR_{1-1/kappa}(eta_L) = inf_{0<=t<=1} t+kappa E[(eta_L-t)_+].
To prove the upper bound, write E[w eta]=t+E[w(eta-t)] and retain the
positive part with w<=kappa. Equality follows by putting weight kappa on
the strict upper tail and using a fractional weight at the boundary atom
to make E w=1. The fractional boundary is essential for discrete laws.
Conditional Jensen gives E[(eta_L-t)_+] <= E[(Xbar-t)_+], with
Xbar=(X1+X2)/2. Thus CVaR(Xbar) is a valid, possibly noisy relaxation.

Evaluate two CVaR confidence procedures, separately at delta=.001:

1. **DKW:** epsilon=sqrt(log(2/delta)/(2M)). On the uniform CDF event,
   E[(Xbar-t)_+] <= empirical_mean[(Xbar-t)_+] + epsilon*(1-t)
   for every t in [0,1]. Minimize t+kappa times this upper expression at
   all empirical breakpoints and endpoints; cap at 1. The single uniform
   event permits minimization without another selection penalty.
2. **Split empirical Bernstein:** first M/2 independent views choose t as
   an empirical quantile at 1-1/kappa. The other half bound the mean of
   (Xbar-t)_+ in [0,1-t] with one-sided empirical Bernstein at delta.
   Use t+kappa times that mean bound, capped at 1. The first half is used
   only for t. This is a separate valid method, not a minimum with DKW.

Two simpler variance comparisons use the same archived Xbar observations:

3. **Mean only:** bound mu with a two-sided empirical Bernstein interval
   [L,U], delta/2 per tail. Bound V<=max_{u in [L,U]} u(1-u), using 1/2
   clipped to that interval. Combine U+sqrt((kappa-1)V), kappa U and 1.
4. **Unpaired second moment:** allocate delta/2 to the two-sided mean and
   delta/2 to an upper mean of (Xbar-b)^2 in [0,max(b^2,(1-b)^2)].
   Subtract distance(b,[L,U])^2, clip to [0,1/4]. Since Var(eta_L) <=
   Var(Xbar), this bounds the desired V without the paired-product identity.
   Use the same joint-event minimum as the paired procedure.

The classical CVaR minimization, including atoms, follows Rockafellar and
Uryasev (2002), whose primary author manuscript Section 3 and Theorem 10 were
read selectively: https://sites.math.washington.edu/~rtr/papers/rtr187-CVaR2.pdf .
The statement of the one-sample DKW--Massart bound is checked in the
introduction of Wei and Dudley's primary paper, https://arxiv.org/pdf/1107.5356 .
Massart's original proof was not accessed, and the latter paper's two-sample
results are not used. None of these risk/concentration tools is claimed new.

## Reporting and checks (V4--V5)

Retain all new bounds and 480 binomial projections, alongside the previously
saved three procedures. Decompose each paired upper variance into signed raw
centered-product mean, square-root concentration term, range term, distance
subtraction and clipping. Negative raw estimates are retained. Report paired
critical-count differences at each fixed test size. Held-out pointwise intervals
measure single-image frequency uncertainty only, not calibration variability.
One calibration batch per stack does not characterize that variability.

Check the CVaR population formula against a separate finite-distribution linear
program, the DKW objective against explicit breakpoint enumeration, and Jensen
ordering on an enumerated conditional-noise model. Add an exact noise-adaptive
amplitude counterexample and near-degenerate cubic checks. The existing arrays
cannot recover optimal replication allocations; no best-allocation or
calibration-repetition result is asserted. Such a study requires fresh draws.
