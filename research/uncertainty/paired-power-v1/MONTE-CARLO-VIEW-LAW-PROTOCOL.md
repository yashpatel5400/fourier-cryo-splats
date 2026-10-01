# Candidate validation under a bounded viewing law

1 October 2026 UTC. A distinct, narrower model follows the failed deterministic
covering-cost gate. No arbitrary-view guarantee is asserted. This is a
feasibility experiment using fixed synthetic cell densities derived from real
reference maps, not experimental density coverage or a new reconstruction.

## Statistical statement and proof

Fix a candidate density, transferred continuous Fourier-slice simulator m(R),
known independent proper Gaussian noise, amplitude interval A, and a real
translation-invariant score f. The score and threshold c must be chosen
independently of calibration and tested images. Let Q be a specified viewing
distribution, here Haar on SO(3). Assume each tested particle is independent
and its viewing law P_i satisfies dP_i/dQ<=kappa, with kappa>=1 known. The
amplitude may vary within A and depend on orientation. Translations can vary
independently of noise: their phase action leaves the score invariant and
preserves the proper Gaussian noise law. This does not cover unknown windows,
correlated or noncircular noise, CTF errors, or an uncalibrated viewing bound.

For an independent simulation draw R~Q and epsilon from the noise law, form

    Z = 1{ max_{a in A} f(a m(R)+epsilon) > c }.

The maximum is over a scalar cubic polynomial, found from its endpoints and
real stationary points. Let K be the sum of M iid draws of Z, and let U(K)
be the one-sided Clopper--Pearson upper probability at confidence 1-delta.
Define p_bound=min(1,kappa U). For n tested particles, reject when the number
with f(Y)>c reaches the smallest integer k satisfying

    Pr{Binomial(n,p_bound)>=k} <= alpha-delta.

The false-rejection probability, jointly over calibration and test data,
is at most alpha, in exact arithmetic under the stated model. To prove this,
the noise-wise amplitude envelope dominates each possible amplitude's event.
Changing the viewing law increases its probability by at most kappa. Except
on an event of probability delta, the simulation probability is bounded by U.
Conditional on that event and the independent calibration, independent test
indicators have probabilities at most p_bound; their sum is stochastically
dominated by the displayed binomial variable (couple independent uniforms).
A union bound gives alpha-delta+delta. Equal particle viewing laws are not
needed; independence and the common probability bound are. This statement
does not condition on a realized arbitrary set of poses.

The result is a composition of classical binomial confidence and stochastic
domination, not an established novel theorem. General Monte Carlo testing
and nuisance maximization long predate this experiment. Dufour (2006),
[author's primary PDF](https://jeanmariedufour.research.mcgill.ca/Dufour_2006_JE_MCT.pdf),
DOI 10.1016/j.jeconom.2005.06.007, was read selectively through its introduction
and Section 2.1: it gives simulation-rank tests and distinguishes finite-sample
maximization from approximate nuisance substitution. We have not replicated
its other proofs. Berger and Boos (1994),
[publisher abstract](https://www.tandfonline.com/doi/abs/10.1080/01621459.1994.10476836),
is an additional prior-art lead on confidence-set nuisance maximization;
only its primary abstract/search excerpt was accessible here. The present
event-envelope calculation is not claimed to be either paper's algorithm.

## Locked simulation

Use all twelve saved removal directions (power / power+bispectrum; amplitude
{1} / [.9,1.1]; three stacks). Do not refit directions or select successful
stacks. Four 10076 zero directions give a deterministically constant statistic
and are retained explicitly as uninformative, without expensive projections.
The scalar threshold is each direction's previously saved noiseless
alternative mean; its finite 4,160-view origin is declared and poses no
calibration reuse because it is fixed before these random draws.

For each nonzero stack, draw 131,072 calibration orientations/noises and
131,072 independent held-out orientations/noises. Each orientation is a fresh
continuous Haar draw, not sampling a saved finite catalog. Generate the true
and region-removed physical-cell projections with the exact cell Fourier
operator (NUFFT evaluated in floating point); use the first saved transfer
profile and unit noise per real/imaginary coordinate. Process batches of 512.
Seeds are 261012+dataset+stage*1,000,000. Store the RNG state, hash input files,
and retain scalar scores for every simulated particle and method. Within a
stage use common rotations/noise across methods and true/removed maps; the
two stages use independent seeds. No interpolation guide evaluates a score.

Calibrate both removed-candidate and true-map control envelopes. On the
held-out batch record true-map amplitude-one scores, removed-map amplitude-one
scores, and both amplitude-envelope scores. The true-map candidate is a
correct-null control evaluated on the true-map held-out data. Use alpha=.05,
delta=.001; report kappa=1,1.01,1.1,2 and n=1,000,10,000,100,000 separately.
These bounds are assumptions, not empirically calibrated preferred-view limits.
There is no unadjusted union of methods, thresholds or candidates.

Report single-particle event frequencies with exact pointwise 95% binomial
intervals. Transform their endpoints through the binomial rejection tail to
give projected n-particle rejection probabilities and pointwise intervals,
conditional on the independent realized calibration. These are power/type-I
projections, not observed repeated-dataset rejection rates or simultaneous
confidence bounds across the table. Include independent one-particle checks
of the null event envelope and all failures. Runtime, numerical operator
checks and random seeds accompany the output. Independently sum all physical
cells for the first orientation of each stage and each nonzero stack, for
both maps, and require agreement within 1e-8 in each transferred Fourier
coordinate. Preserve those orientations and means. Ordinary floating-point
evaluation remains distinct from certified numerical error control.

No rental or experimental-data claim follows from this study. Advancing
requires useful performance, sensitivity to realistic preferred viewing,
known simulator error/noise and a credible contribution beyond classical
simulation testing. A positive Haar-only outcome does not repair the existing
three full-review rejections.
