# Shift-invariant third-moment feasibility gate

1 October 2026 UTC; post-outcome method development. Earlier quadratic
directions have verified continuous-pose counterexamples. This gate asks
whether retaining translation-invariant phase information helps distinguish
the same local-removal candidates. Bispectra, moment matching and Gaussian
Hermite variance calculations are established tools, not claimed inventions.
Relevant prior art is recorded in `../MOMENTS-AND-IDENTIFIABILITY-READING.md`;
Perry et al., SIAM J. Math. Data Sci. 1(3):497–517 (2019), further distinguishes
individual alignment from collective signal recovery. A useful new uncertainty
procedure is still missing.

## Statistics and exact conditional noise calculation

For unique nonzero Fourier coordinates with no opposite pairs, choose 512
distinct triples a<b with qa+qb=qc and c distinct from a,b, without replacement
using seed 261010. Selection uses frequencies only. For Y=m+epsilon with
independent N(0,v) real and imaginary noise coordinates, use

    (|Yj|²-2v)/2, Re(Ya Yb conjugate(Yc))/2,
    Im(Ya Yb conjugate(Yc))/2.

Their expectations are the corresponding noise-free features. Continuous
translations cancel exactly in every product. This presumes the stated Fourier
sampling, noise and CTF model; windowing or correlated/noncircular noise can
break these assumptions. Unit pure-noise marginal variances do not imply
independent coordinates at signal or a fully whitened statistic.

For a fixed real contrast f of these coordinates, the exact conditional
Gaussian variance is

    v ||grad f(m)||² + v² ||Hess f(m)||F²/2
       + v³ ||D³ f||F²/6.

Proof: the polynomial has degree at most three. Its Laplacian is constant,
so its linear Gaussian-Hermite coefficient is grad f(m); its quadratic and
cubic coefficients are the ordinary derivative tensors divided by 2 and 6.
Orthogonality of different Hermite orders gives the identity. All overlapping
triad terms are summed before squaring. This is the classical finite Hermite
expansion specialized to this statistic, not a new general theorem. Independent
six-dimensional four-node Gauss-Hermite integration checks mean and variance;
a separate simulation checks shared-frequency cross terms.

## Locked oracle experiment

Use all 4,160 saved orientations, the first original CTF/noise profile, and
the true and fully removed maps on 10028/10049/10076. The alternative is the
true map's uniform empirical distribution on those orientations, not an exact
Haar integral. Compare power alone with power plus 512 complex bispectra.
For each, fit the candidate moment hull at amplitudes {1}, and {.9,1,1.1}.
This gives 12 removal fits and 12 exact true-map controls. Amplitude ranges
are sensitivity assumptions, not experimentally calibrated bounds.

The optimizer projects the alternative moment vector onto the convex hull
of null moment vectors. A feasible simplex mixture supplies an upper distance.
For unit h, h'target-max_l h'atom_l supplies a lower distance. Preserve both,
their gap, all mixture weights and coefficients. Accelerated projected gradient
uses checked backtracking, at most 2,500 steps and 90 seconds per fit, with
distance-gap tolerance 1e-5. A limit is not convergence or nonmembership.

After fitting, profile the contrast mean and conditional variance analytically
over the full corresponding amplitude interval: both are polynomials of degree
at most four, so evaluate endpoints and all real stationary points. Check the
10,000 saved independent orientations from the full-frequency gate. These
orientations have been used in previous development but are not optimizer atoms
for this gate. Save each maximum and changed separation. They are still a finite
set, not a certificate on SO(3).

Report the alternative's exact noise variance averaged across its 4,160 means
plus the variance of its conditional means. For the enlarged finite catalog,
report the largest conditional noise variance over rotations and continuous
amplitudes. A conservative finite-catalog sample-size diagnostic follows from
Cantelli: for separation delta>0, independent particles, alpha=.05 and beta=.2,

    n >= (sqrt(19 Vnull) + 2 sqrt(Valternative))² / delta².

To verify: condition on every null pose/amplitude. Each conditional mean is at
most the null maximum and each noise variance at most Vnull. Cantelli controls
the sample-mean test at alpha with threshold nullmax+sqrt(19 Vnull/n).
Under the specified iid alternative, a second Cantelli bound gives miss rate
at most beta. This is a sufficient sample size under the stated finite-catalog
model, not the necessary sample complexity or measured rejection power. No
unbounded-amplitude, unknown-covariance, continuous-view or density-coverage
claim follows. If the checked separation is nonpositive, record no such bound.

Do not exponentiate a cubic Gaussian contrast without a valid normalizer:
its moment-generating function need not exist. This gate is not the earlier
bilinear e-value construction. No experimental particles or new noise trials
are used to fit the contrasts. Source/protocol hashes are committed before the
first molecular outcome, and all failures/limits remain visible.
