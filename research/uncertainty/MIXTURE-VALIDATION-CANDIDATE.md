# Candidate: continuous-pose structural validation by likelihood envelopes

30 September 2026. Prospective development, separate from the manuscript's
completed conditional linear-feature studies. No empirical accuracy, power or
novel general inference theorem is claimed. The current cubic fits and RELION
comparisons continue unchanged. This candidate addresses their unresolved
dependence on externally supplied small pose balls; it changes the estimand
to structural compatibility, not arbitrary local-density confidence.

## Model and established statistical foundation

For a candidate density theta, independent validation images have Gaussian
densities `p_i(y | theta, R, nuisance)` conditional on their known imaging
operators. The viewing distribution G is unknown but shared across images;
its independence from imaging parameters is an assumption. Image-specific
translations may instead be profiled over an explicitly declared set. A
heterogeneous extension requires specifying the ensemble in the null, not
calling pooled EMPIAR-10076 homogeneous.

An independently fitted, normalized predictive density q on the whole
validation sample gives the classical split-likelihood e-variable
`q(Y) / sup_G,nu p_theta,G,nu(Y)`. Any pointwise upper bound on the denominator
remains valid, including one computed using Y. Markov's inequality gives
type-I error at most alpha for rejecting when the ratio is at least 1/alpha.
Inverting the tests covers the true candidate when it is in the declared
family. A finite list is not a confidence set for every continuous density.
The numerator must be a normalized density fixed independently of validation
noise; an unnormalized profile score is insufficient. This is universal
inference (Wasserman, Ramdas and Balakrishnan, 2020), not a new principle.

## A finite convex upper bound for the continuous viewing law

Partition the whole allowed rotation domain into cells C_c. Suppose a
computable envelope satisfies

`U_ic(y_i) >= sup_{R in C_c, allowed nuisance} p_i(y_i | theta,R,nuisance)`.

The true cell masses w_c=G(C_c) lie on a simplex, so the joint log likelihood
is at most `max_w sum_i log(sum_c U_ic w_c)`. This is a concave maximization
with linear simplex constraints. It is an upper relaxation of the continuous
model: the maximum inside a cell need not occur at the same rotation for
every image. Refining cells can reduce this relaxation.

For any positive anchor z_i, the tangent inequality for log gives the upper
bound

`D(z) = sum_i log(z_i) - n + max_c sum_i U_ic / z_i`.

Indeed, `log x <= log z + x/z - 1`; maximize the resulting linear function
over the simplex. With z=Uw at a feasible w, the primal-to-dual gap is
`max_c sum_i U_ic/z_i - n`. This is a standard mixture-likelihood dual
certificate. An EM weight iteration can supply anchors; convergence is not
required for a valid real-arithmetic upper bound. A local density or pose
optimization alone cannot replace U or D.

## A Gaussian cell envelope and nuisance variance

For d real whitened coordinates, let every mean in a cell lie in a Euclidean
ball of radius b around m. Reverse triangle inequality gives residual lower
bound `r=max(||y-m||-b,0)`. For sigma in a declared positive interval,

`log U = max_sigma [-d log(sigma sqrt(2pi)) - r²/(2 sigma²)]`.

The maximizing sigma is `r/sqrt(d)` clipped to the interval. If the interval
includes arbitrarily small sigma and r=0, the envelope is infinite; the test
must retain that uninformative outcome. Profiling variance separately within
each envelope is valid even for a shared true variance, but can be loose.
It does not cover unknown colored covariance merely by renaming a scale.

Fourier Gaussian splats permit analytic mean evaluation and derivatives.
A future continuous rotation-cell implementation must bound their image-space
variation, cover the entire allowed rotation domain, include translations and
CTF uncertainty as claimed, and retain integration/floating-point limitations.
Those pieces are not yet implemented. A discrete orientation grid with zero
cell radii is a different, restricted model, not a continuous-pose guarantee.

## Falsifiable development gate

First verify the algebra against an independent finite conic mixture fit and
direct Gaussian envelope calculations, including degenerate/infinite cases.
Then declare a bounded feasibility study before generating outcomes. Start
with an explicitly discrete orientation model and an optimistic independently
specified predictive density. If even this test cannot separate meaningful
structural changes at relevant noise levels, do not expand to an expensive
continuous-pose implementation. If it can, the next required measurement is
the loss of power and runtime caused by certified cell envelopes, unknown
noise scale and realistic acquisition variation. Discrete-model success alone
must not enter the paper as a continuous or experimental validation result.

Closest prior art includes BioEM, CryoLike, moment-based structural metrics,
nonuniform/SubspaceMoM, and universal inference; see the targeted reading
notes. A publishable contribution would require more than the elementary
inequalities above: useful continuous-pose computation and persuasive
experimental structural validation remain unestablished.

## Development status after the initial proposal

The restricted discrete study and shared common-scale follow-up are complete.
Continuous SO(3) enclosures, scaled classical mixture bounds and a local
curvature refinement are now implemented and independently checked. Thus the
future-tense implementation paragraph above records the proposal at its start,
not today's status. All three initial continuous calculations exhaust their
split budgets with very large likelihood brackets. Envelope-weight optimization
reduces the gaps but leaves 17,077--20,253 log units above the best recorded
feasible values. The results note preserves these negative outcomes. A new
full-cover refinement protocol is running on the same observations. Unknown
translations, continuous shared-scale optimization, a practical independently
learned normalized numerator and experimental calibration remain unfinished.
