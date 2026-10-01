# Descriptive background-spectrum inventory

**Post-outcome metadata correction:** the following frozen protocol includes an
incorrect inference about source identities. The existing metadata joins recover
229 / 137 / 351 source groups. See the [dated correction](BACKGROUND-SOURCE-GROUP-CORRECTION.md).
The numerical protocol and its original wording below are preserved.

1 October 2026 UTC. Use all 8,192 previously downloaded 64×64 particle images
per stack, without additional selection. Analyze all four nonoverlapping 8×8
corner patches, pooled and separately, in the full cohort and both recorded
source halves. Report second moments, covariance and fourth moments, not a
pure-noise confidence interval. Molecular signal, neighboring particles, ice,
source preprocessing and shared micrographs can contaminate these patches.
The pose/CTF metadata does not identify all source micrographs: the selected
blob paths number 229 / 2 / 1. Do not treat those as independent acquisitions.

The downloader Fourier-crops to 64 pixels and zeros the output Nyquist row
and column; otherwise it neither whitens nor normalizes the images. A valid
white-noise reference must include that operation. For periodic coordinate
lag h, the resulting one-dimensional normalized covariance is
(64 1{h=0 mod 64} − (−1)^h)/63. Its tensor product gives the 8×8 patch
covariance. This identity can be checked against the explicit retained Fourier
basis, avoiding a spurious report that removed Nyquist modes are noise color.

Subtract each patch's mean, transform by the orthonormal 2D DCT-II, remove the
DC coordinate, and whiten only by this analytically specified *reference*
covariance using Cholesky. The operation removes known preprocessing effects;
it does not empirically whiten the real data. Save the 63-dimensional second
moment and centered covariance, normalized by mean second moment, per group.
Retain diagonal extrema, spectral eigenvalue extrema, off-diagonal Frobenius
fraction, mean-vector energy, kurtosis and particle-energy quantiles. Save the
unwhitened DCT power and reference diagonal as well. Record the full-cohort
and half-specific pooled power ratios without interpreting halves as independent
noise exposures. No threshold or test is selected from this inventory.

Verify the reference covariance by an explicit complex Fourier basis, and the
reported selected patches by an independent cosine-transform matrix, centering
and solve. This gives a descriptive answer to whether the archived pixels
look compatible with precisely known white noise. It cannot measure pure
noise covariance inside the molecule or certify sub-percent noise knowledge.
