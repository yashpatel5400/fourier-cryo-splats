# Frozen additional-exposure prediction protocol v1

Date: 29 September 2026. This protocol is frozen before downloading or evaluating
its selected particle pixels. Its purpose is a conditional experimental-image
prediction check, not calibration of density intervals, new-structure
extrapolation, or independent-pose gold-standard validation. The existing
8,192-image pool per dataset is development data and remains labeled that way.

## Fixed cohort

Use EMPIAR-10028, 10049 and 10076. In each dataset, exclude **every source
exposure group** represented in the entire development pool. Group identities
are joined using source image identities for 10028/10049 and the deposited
Frealign FILM field for 10076. Apply the published cryoDRGN particle filters and
previously documented missing-source exclusions. Permute contiguous source
blocks of 64 with seed 609421 + integer EMPIAR accession; retain the first 4096
eligible particles, then sort indices. Metadata-only selection is already
recorded in the CSVs and manifests: 115, 71 and 185 new groups respectively,
with zero development-group overlap. No selected image pixel has been inspected
for this protocol at freezing.

Download exactly the recorded source indices into separate data directories,
using verified HTTP byte ranges and the existing Fourier crop to 64 pixels.
Record headers, payload SHA256s, source dimensions and output-array hashes.
Retries may replace failed transfers of the same indices. If a source becomes
unavailable, report incompleteness; do not select replacements based on results.
Consensus poses/CTFs were published using larger datasets, so they remain
conditional inputs even though these images are new to our fitted models.

## Locked models and evaluation

Use the already trained Gaussian and voxel ridge reconstructions, averaged over
the two development exposure-group halves. Use the two stock cryoDRGN homogeneous
networks averaged in prediction space, at epochs 60, 50 and 100 for 10028, 10049
and 10076. Those epochs were selected only by development tuning-group NMSE.
The frozen model manifest records coefficient/checkpoint/config hashes and the
existing classical training scale. No model, target, normalization scale or
epoch is reselected using the additional exposures. Per-image background
normalization and the fixed spherical window are part of the prespecified
observation transformation, identical across evaluated methods.

Evaluate all 1410 conjugate-pair representatives inside radius 30 of the 64
pixel grid, with the same published translation centering and CTF. Training
objectives differ (stock neural radius 32 versus classical radius 30); retain
this limitation. Check neural direct predictions against independently saved
3D-Fourier map values using the established five-point check.

Primary result: total squared complex Fourier prediction error divided by total
observed Fourier power (NMSE), for all three methods in all three datasets.
Also report normalized prediction/observation correlation. Retain per-particle
errors and power locally and per-exposure aggregates publicly.

Use 10,000 paired whole-exposure bootstrap resamples, with seed 609422 + integer
accession. Report Gaussian-minus-neural and voxel-minus-neural NMSE differences,
95% marginal percentile intervals and 99.1667% Bonferroni-adjusted percentile
intervals for the six prespecified contrasts. These are ordinary group-bootstrap
approximations conditional on the supplied poses, not finite-sample coverage
proofs. No cross-dataset pooling, best-dataset selection or stopping for a
favorable difference. Reconstruction FSCs do not change because models are
frozen; do not imply that this new-particle experiment creates new half maps.

## Integrity and deviations

Commit this protocol, the selection and model hashes before the first download.
Archive the exact executed source bytes before evaluation. Verify both model
and cohort hashes in the evaluator; fail on drift. Any necessary numerical bug
fix is a dated protocol amendment with the original output/error retained;
post-outcome scientific changes require a separately labeled new study.
Run all three datasets irrespective of outcome. The results provide one frozen
confirmation component of the larger ongoing study, not its completion.
