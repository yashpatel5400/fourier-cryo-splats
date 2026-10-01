# Simulation-only planted-truth audit

1 October 2026. This is a post hoc diagnostic of saved importance samples,
not a modification of either frozen numerical gate and not a new estimator
for experimental images.

## Elementary identity and scope

Fix y and an integrable positive kernel f_y(x), with integral I_y. Let q_y
be a normalized proposal with positive density wherever f_y is positive.
Conditionally on y, let X_1 have density f_y/I_y and X_2,...,X_n be independent
q_y draws, also independent of X_1. Define w(x)=f_y(x)/q_y(x) and
I_star=(sum_i w(X_i))/n. Then

    E[1/I_star | y]
      = (1/I_y) E_(q_y^n)[n w(X_1)/sum_i w(X_i)]
      = 1/I_y.

The last equality follows by exchangeability: summing the n equal expectations
gives E[sum_i w_i/sum_i w_i]=1. The same argument works for a fixed replaced
index because the statistic is symmetric. It does not require a random index.
The other n-1 points must remain fresh draws from the full proposal. Merely
appending a truth point without accounting for n, or replacing a selected
large/small weight, would be a different calculation.

For a correctly generated joint sample (X,Y), X conditional on Y is a posterior
draw. The primary source [Grosse, Ghahramani and Adams (2015), Section 4.2](https://arxiv.org/html/1511.02543v1#S4.SS2)
uses this generative/posterior identity for reverse Monte Carlo likelihood
diagnostics. We read Sections 3.1 and 4.1–4.2; we do not implement that paper's
AIS/SMC algorithms. The finite-i.i.d.-importance argument above is stated
and proved directly, not claimed as a new result.

The proposal may depend on y and an independent fixed training catalogue, but
must not be fit using the generating pose beyond information in y. The exact
posterior statement averages over the matched model's generative joint law;
it does not supply repeatable conditional coverage for an already fixed true
pose. In this experiment only the actually generating state has that premise.
It fails for the alternative map on the same image and for experimental data.

For independent simulated images, the product of forward importance estimates
is unbiased for the stack integral, and the product of reciprocal planted
estimates is unbiased for its reciprocal. Markov and a union bound give

    sum log I_forward - log(1/delta_L) <= sum log I
        <= sum log I_star + log(1/delta_U)

with joint probability at least 1-delta_L-delta_U. This is a stochastic bracket,
not a deterministic numerical certificate. We report delta_L=delta_U=.025
separately for each fixed stack/state/bank. There is no simultaneous claim
across reported brackets. Paired null/altered states share noise and poses and
are **never** multiplied together. The two banks share the same planted pose
and are not treated as independent replicates of an upper bound.

## Fixed retrospective calculation

For all 128 saved image cases in each stack and both 8,192-point banks:
replace draw zero, and only draw zero, with the saved generating rotation.
Use the matching generating state and physical kernel -||noise||^2/2.
Reconstruct every saved observation from archived generating means/noise.
Retain original and planted log integrals and their signed difference.
For the first proposal use its full normalized density; for the catalogue
repair use its own full density at truth with all catalogue centers.

Also retain matching-state weight mass on Haar-labelled points, whether its
maximum weight comes from Haar, both integral delta-method SEs, the shared-draw
log-ratio delta-method SE, and the between-bank difference divided by its
estimated SE. These SEs are plug-in diagnostics, unreliable in poorly sampled
tails; they are not certified error bounds. Curvature clipping is determined
from all original raw Hessian eigenvalues and the original .25–30 degree
cutoffs, with original mode weights, without refitting.

Report every image and stratify only by the original predeclared .01 ratio
criterion. This is post-outcome failure diagnosis, not predictive validation.
It cannot tune the running catalogue proposal or permit a third revision.
