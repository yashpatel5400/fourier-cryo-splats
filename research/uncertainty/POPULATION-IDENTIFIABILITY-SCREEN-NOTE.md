# A population target needs an identifiability screen first

1 October 2026 UTC. Development mathematics prepared while the focused methods
consultation runs. This is not a selected final method, a new statistical
theorem claim or evidence that the real-data assumptions hold. It spells out a
possible falsification check before implementing population uncertainty.

## Target and elementary support characterization

Fix two templates rho_0,rho_1 and one known imaging operator family. Let
m_s(xi) in R^d be the noiseless observed coordinates for state s and nuisance
xi (including orientation and any allowed amplitude/shift). Let M_s be the
compact image of a compact nuisance space under a continuous forward map.
Allow arbitrary state-dependent probability measures on that nuisance space.
The target pi is the fraction in state 1 among the sampled particles, not an
equilibrium free energy or the pre-selection biological population.

Assume Y = M + epsilon with independent, known, nonsingular Gaussian noise.
Write mu for the probability law of the noiseless observed mean M. In the
infinite-data population model, convolution identifies mu uniquely: the
characteristic function of a Gaussian is nonzero everywhere, so divide the
observed characteristic function by it and use uniqueness of characteristic
functions. This does not identify individual orientations and is generally an
ill-conditioned operation with finite data.

If mu is supported on M_0 union M_1, the compatible state-1 fractions are

    mu(M_1 minus M_0) <= pi <= mu(M_1).

For finite nuisance catalogues this characterization is exact and elementary.
For continuous spaces it holds assuming measurable probability lifts from
each mean support to its allowed nuisance space; this assumption is explicit
here rather than relying on an unproved selection result.

Proof: mass supported only on M_1 must belong to state 1, and no state-1 mass
can lie outside M_1. Mass on the intersection can be assigned arbitrarily
between states. Splitting that intersection mass in the desired fraction and
lifting each resulting measure supplies a compatible joint state/nuisance
law. Independent Gaussian convolution leaves the observed law unchanged.
Thus disjoint image supports identify pi in this population model. An exact
shared mean permits an observed distribution compatible with every pi in
[0,1]. Priors, density floors, state-view independence or smoothness constraints
can change the identified set; none may be silently added or removed.

This statement is a classical Gaussian-deconvolution/support argument, not a
claim that population inference with arbitrary preferred views is novel or
impossible for every molecule. An overlap of two finite catalogues is not the
same as an overlap of the continuous image supports. Failure to find an
overlap also does not certify their separation.

## A verifiable finite-sample obstruction

Choose any admissible nuisance xi_0,xi_1, put a=m_0(xi_0), b=m_1(xi_1), and
let delta=||Sigma^(-1/2)(a-b)||. Consider two experiments: every particle is in
state 0 at xi_0, or every particle is in state 1 at xi_1, with n independent
Gaussian observations. Their KL divergence is n delta^2/2. Pinsker's inequality
therefore bounds total variation by min(1, sqrt(n)*delta/2).

If a random confidence interval C has coverage at least 1-alpha in both
experiments, then under the first experiment

    P({0,1} subset C) >= max(0, 1-2 alpha-min(1,sqrt(n)*delta/2)).

Proof: under the first experiment, inclusion of 0 has probability at least
1-alpha. Inclusion of 1 has probability at least 1-alpha-TV by transferring
the second experiment's coverage event. The union bound gives the displayed
inequality. An interval containing both endpoints has diameter at least one.

Any numerically found admissible pair supplies an upper bound on the distance
needed in this obstruction, after independently verifying its physical means
and numerical error. A global minimum is not needed for such a negative
witness. Conversely, a failed local search supplies no positive identifiability
certificate. The bound can be zero and therefore uninformative; report that
outcome without calling it an impossibility result.

The point-concentrated viewing laws are allowed only in the arbitrary-view
model. If the intended model enforces a positive uniform component or a
restricted acquisition distribution, these witnesses must be replaced with
admissible distributions before invoking the bound. Approximate physical means
must include their errors in delta. Shared preprocessing, dependent images,+unknown colored noise and imperfect templates invalidate the simple
experiment unless separately modeled.

## How this would constrain a new method

The next population proposal would need to distinguish structural ambiguity,
finite-sample uncertainty and numerical approximation. Profiling joint state
and view weights over a fixed catalogue is a classical convex-mixture problem;
its interval cannot automatically be treated as continuous-pose coverage.
An externally validated template pair and a measured nuisance model are needed
before running a scientific usefulness study. The CAHRA high-signal convention
supplement is suitable for coordinate checks, not for establishing realistic
population calibration or resolving arbitrary-view identifiability.

No numerical population screen has been run for this note. No favorable
separation, confidence interval, new theorem or population result is claimed.
