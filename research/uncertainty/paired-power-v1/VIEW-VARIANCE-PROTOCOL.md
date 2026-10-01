# Conditional-noise grouping and viewing-distribution sensitivity

1 October 2026 UTC. Post-outcome development after v0.7.3-dev. That release's
scalar kappa-times-probability bound is ineffective on 10049 at kappa=1.1.
This does not show that preferred views actually destroy discrimination.
The new gate uses repeated conditional noise to separate variation between
views from noise variability, while retaining the same uncalibrated noise,
CTF, amplitude and viewing-density assumptions. No novelty or experimental
coverage is claimed from composing the classical inequalities below.

## Amplitude envelopes with repeated simulated noise

Partition [.9,1.1] into 16 closed equal cells. At one view R, draw L=32
independent noises and define Z_lj as the indicator that the cubic contrast
exceeds the frozen threshold somewhere in cell j. Both endpoints and real
stationary points give the exact polynomial cell maximum in real arithmetic.
Set X=max_j (1/L) sum_l Z_lj. The old per-noise envelope would instead be
T=(1/L) sum_l max_j Z_lj, and X<=T for every simulation.

For any fixed physical amplitude a, the cell containing a supplies
1/L sum_l 1{f(a m(R)+epsilon_l)>c}<=X pointwise. Consequently
eta(R)=E_noise[X|R] dominates the tail probability at every a. This holds
even if amplitudes vary across actual particles and depend on orientation.
The grouping is only a simulation construction, not an assumption that
experimental particles share an amplitude. E max differs from max E; the
finite-noise maximum retains an upward conservatism that is disclosed.

Draw two independent noise groups at each R. For a center b fixed from the
older independent calibration, E[(X1-b)(X2-b)|R]=(eta(R)-b)^2. This removes
the within-view noise variance from the second moment. Let mu=E_Q eta and
V=Var_Q eta. A density ratio w=dP/dQ bounded by kappa has E_Q w=1 and
E_Q(w-1)^2<=kappa-1. Cauchy--Schwarz then gives

    E_P eta <= mu + sqrt((kappa-1)V).

This bound concerns the conditional mean eta, not the variance of a single
noisy indicator or of X. The simpler kappa*mu bound remains valid too.

## Finite Monte Carlo confidence allocation

The independent sampling unit is a view with its two noise groups, not an
individual noise replica. Use an empirical Bernstein radius for iid
Z_i in [l,u]: sqrt(2 s² log(2/d)/M)+7(u-l)log(2/d)/(3(M-1)), with unbiased
sample variance. The primary COLT 2009 paper by Maurer and Pontil,
[Theorem 4 and Section 2](https://www.learningtheory.org/colt2009/papers/012.pdf),
was read selectively, including Theorem 11's derivation from the variance
bound and Bennett's inequality. Its learning-theory experiments and other
proofs were not replicated. We use this established bound, not claim it as
a new concentration result. Rescaling handles the signed product below.

Allocate delta/2 to a two-sided interval [L_mu,U_mu] for the mean of
(X1+X2)/2 (delta/4 per tail), and delta/2 to an upper bound U_Y on the mean
of Y=(X1-b)(X2-b). Its range is [-b(1-b),max(b²,(1-b)²)]. Then, except
with probability delta,

    V <= V_U = min(1/4, max(0,U_Y-distance(b,[L_mu,U_mu])²)).

Use p_bound=min(1,kappa U_mu,U_mu+sqrt((kappa-1)V_U)). Both expressions
hold on the same joint confidence event, so this minimum needs no additional
selection correction. The independent-particle binomial domination proof in
the preceding protocol then gives total false-rejection probability at most
alpha using level alpha-delta. Ordinary floating-point evaluation remains
separate from certified numerical error control.

## Declared comparison and computations

Study all four nonzero ranged-amplitude directions: power and power+bispectrum
on both 10028 and 10049. These are selected follow-up cases; four fixed-amplitude
and all four zero 10076 directions stay in the previous complete inventory,
without new fits or claims of a three-stack positive result. The score and
threshold are unchanged. Each candidate's b is its old stage-0 envelope event
frequency, fixed before new simulations. Both true-map and region-removed
candidates are calibrated; true and removed amplitude-one held-out data are
both retained as null/alternative controls.

For each stack generate M=32,768 fresh continuous Haar calibration views,
two independent groups of 32 noises per view, and 131,072 fresh independent
held-out Haar views with one noise each. Use the identical physical cell
densities and first transfer profile. Calibration batches have 32 views;
held-out batches have 512. Seeds are 261013+dataset+stage*1,000,000. Stages
are independent; maps and methods share draws within a stage. Preserve grouped
X1/X2 and T1/T2 arrays, all held-out scalar scores, RNG states and input hashes.
Verify the first view of each map/stage by an independent physical-cell sum.

Compare three prespecified methods at the same simulation budget: the
individual-noise envelope with kappa-times-mean, the grouped envelope with
kappa-times-mean, and the grouped viewing-variance bound. Each direct method
uses its own one-sided mean confidence at delta=.001; the variance method
uses the joint split above. Their guarantees are separate; there is no
unadjusted minimum across methods. Use alpha=.05, kappa=1,1.01,1.1,2,5 and
n=1,000,10,000,100,000. Report every bound and projected rejection probability
with the same pointwise binomial intervals and limitations as before.

The held-out alternative remains Haar; changing kappa changes the protected
null model and critical value, not the simulated alternative. All negative
outcomes remain visible. More realistic viewing/CTF/noise models and a useful
non-oracle candidate construction remain prerequisites for an experimental
claim, regardless of this gate's results. No GPU rental follows from it.
