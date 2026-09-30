# Discrete-view likelihood feasibility screen

30 September 2026, specified before generating this screen's outcomes.
This is a deliberately optimistic development gate for the separately proposed
likelihood-envelope candidate, not a continuous-pose or experimental guarantee.

Use all three existing acquisition geometries and their fixed first-region
study: 128 CTFs, radius-twelve frequencies, supplied Gaussian noise, and the
unit continuous-L2 64-cell deposited-map generator. The **only** allowed
orientations in this screen are the first 64 of those archived orientations.
Generate each of 128 particles independently from a uniform distribution over
that discrete catalog, with independent known white noise. The source CTF
assignment is retained; latent catalog choices are resampled independently of
CTF. Translations are exactly zero. This does not cover SO(3) or infer original
experimental alignments.

Use an oracle normalized numerator: the exact true discrete mixture with its
known uniform weights and noise. It is fixed before simulation. This is an
optimistic feasibility control, not a learned predictor. For each dataset test
four specified candidate maps: the true map; 50% and 100% attenuation of a
Gaussian-shaped region centered at the already locked first pilot region with
sigma 20 Å; and a zero-signal positive control. Multiply cell coefficients by
`1 - fraction * exp(-||x-center||²/(2 sigma²))`; do not renormalize the altered
maps. Record the actual density and forward-signal changes, which may differ
between stacks. Zero-signal detection alone does not justify proceeding.

For each candidate compute three likelihood ratios using the same numerator:
the oracle known-uniform-law denominator; separate best catalog orientation
per image; and the convex upper bound profiling a common unknown viewing law.
The oracle-law comparator assumes that law; the other two permit an unknown
law, with the mixture bound relying on the shared-law model. Retain two noise
settings: the supplied known variance, and a separate unconstrained positive
variance profiled in every image/cell envelope. The latter is a relaxation,
not a fit with one common variance. It can lose power and may be infinite.
Mean-ball radii are exactly zero; no continuous-cell certification is implied.

Use 16 independently generated repeats per stack, seed 810001 plus the numeric
accession. Retain all repeats, latent indices, generated images, predictions,+likelihoods, weights, optimizer gaps and failures. Optimize for at most 2,000
EM iterations with an absolute log-likelihood gap target 1e-4; any unfinished
fit still uses its available upper bound. A 20-minute overall wall budget is
checked between repeats; retain an incomplete record if it expires.
Report log-e-value distributions and counts exceeding log(20), not an asserted
precise power estimate. The matched true-map oracle-numerator case has ratio
at most one pointwise; it is an algebra control, not empirical calibration.

This screen is a stop/go diagnostic, not a method-selection confirmation.
If local attenuation is undetectable even here, do not claim that the more
difficult continuous/noise-unknown problem is solved. If promising, the next
step must explicitly measure continuous-cell relaxation, practical predictive
fitting, computation and acquisition-model sensitivity before a manuscript
contribution is claimed.
