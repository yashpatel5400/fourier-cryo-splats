# Known-map pose information diagnostic

30 September 2026, declared before these outcomes. The uncertainty experiments
assume local pose radii, while the archived metadata supplies no calibrated
radii. This small diagnostic asks how informative their *known-noise simulation
model* is about a single particle's pose, even if its density were known.
It is not experimental pose calibration or another density-coverage test.

For all three stacks use exactly the 128 particles, radius-twelve frequencies,
CTFs and Gaussian noise convention of the locked first pilot-selected feature.
Use two fixed, separately labeled signals: the unit-norm 64-cell deposited-map
generator, and the unit-norm 24-cell Gaussian-pilot representation already used
in those experiments. Do not select between them. Treat each signal as known;
there is no density fitting or favorable pose/noise selection in this diagnostic.

Compute the exact constant-cell Fourier derivative using analytic first spatial
moments, including each cell's extent. Pose coordinates are three degrees of
rotation and two Angstroms of detector translation. Whiten with the supplied
Gaussian standard deviation. For each particle J, store the full local Gaussian
information I=J'J, its spectrum, and (when nonsingular) I^{-1}. Report rotational
and translational RMS scales from the corresponding covariance blocks, plus
largest directional SDs and all coordinate SDs. Also report the rotational
block inverse if translations were known. Label any numerical rank deficiency;
do not substitute a finite pseudoinverse for unidentifiable directions.

This is the standard Fisher-information/Cramer-Rao calculation. Its inverse
is a local variance lower bound under the usual regularity and local unbiasedness
conditions, or the covariance of the linear Gaussian tangent experiment.
It is not a high-probability error radius for a nonlinear, biased, globally
ambiguous or jointly reconstructed estimator. A Bayesian prior can change
posterior variance. These numbers cannot be promoted to uniform pose balls.
The fixed known density makes this an optimistic reference calculation, not
the actual joint density/alignment Hessian analyzed by Rangan et al. (2024):
https://arxiv.org/html/2411.13263v2.

Verify all five derivatives against independent nonlinear constant-cell forward
differences on a small non-axis-aligned geometry. Record source hashes, exact
input identities, noise, units, rank threshold (1e-10 relative to the largest
information eigenvalue in the stated coordinates), every particle and summary
quantiles. The diagnostic does not estimate a noise law from these particles.
