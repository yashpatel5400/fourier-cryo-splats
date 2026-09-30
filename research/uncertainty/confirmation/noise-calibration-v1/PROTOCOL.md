# Fresh exposure noise-calibration check

This protocol reserves previously unanalysed calibration particles before their
pixels are downloaded. It follows the pilot-selected feature development study;
it does not turn its old inference observations into new confirmation data.
Its narrowly defined purpose is to remove reuse of calibration outcomes for the
already declared twelve estimators. It does not validate nuisance assumptions.

## Cohort and locking

Use EMPIAR-10028, 10049, and 10076. Exclude every exposure represented in either
the original 8,192-particle development pool or the additional 4,096-particle
prediction-v1 cohort. Apply the published particle filter and documented source
exclusions. Source metadata leave respectively 736, 141, and 1,530 unused eligible
exposures. With seed 630929 + accession, permute the sorted unique group labels,
take 128, and choose one eligible source particle uniformly within each selected
group using the same RNG. Sort source indices for download. Preserve selection
CSVs and hash manifests. Selection reads source metadata, not new image pixels.

Before any new pixels are downloaded, freeze all twelve completed fixed-weight
arrays from pilot-selected-fixed-v1 and their target/pilot/amplitude inputs.
Commit their checksum manifest, this protocol, source selections, and the
complete evaluation code. Exact arrays join the next versioned data release.
The same fixed weights are used at every pose budget. Pose certificates may
finish afterwards because their declared algorithms use geometry and fresh
independent random seeds, never calibration pixels. Retain all earlier
experimental outcomes based on reused calibration; this is an additional check.
Do not choose weights, stop fitting, select features or revise budgets in response
to this new calibration. No cohort replacements if a transfer fails: retry the
same verified HTTP byte ranges or report an incomplete dataset.

Download into data/uncertainty/confirmation/noise-calibration-v1/DATASET using
the existing 64-pixel Fourier crop with Nyquist handling. Record the source
range/header checks and exact cached-array hashes. These are extracted particle
images, not detector movies. Published full-data poses/CTFs remain conditional
inputs even though these calibration exposures were unused by our fits.

## Fixed analysis and reporting

Use all 128 calibration representatives per stack. Read uncentered Fourier
coordinates at radius 12; divide by the unchanged independent old-pilot
amplitude. Pull each locked estimator's weights into the same raw coordinate
frame. The original 128 inference particles and exposure groups are unchanged.
Use the grouped variance upper bound, multiplying inference weight row i by
sqrt(n_g(i)); this allows arbitrary jointly Gaussian within-exposure dependence.
It still assumes independent exposures, one common marginal raw-coordinate
covariance across inference and calibration, and weights/design independent
of the relevant noise. Neither metadata disjointness nor this check proves those
conditions. Arbitrary calibration means are allowed by the trace lemma.

For every feature use beta_calibration=0.005/12, delta_spectral=0.000001/12,
and alpha_noise=(0.045-0.000001)/12. Report the same four pose classes: fixed poses;
0.5 A shifts only; 1 degree plus 0.5 A; 2 degrees plus 0.5 A. Reuse the already
declared continuous bias bounds and no-data fallback. The family has twelve
features per pose class. Do not pool pose classes or select a favorable one after
seeing outcomes. Report all 48 intervals, variance bounds, fallback decisions,
zero exclusion and approximate deposited-map inclusion. Report the ratio of
fresh to old calibration SD bounds for every feature, irrespective of direction.
No map-inclusion statistic is a density-coverage estimate.

A successful computation is not a successful validation of the scientific
assumptions. Uniform pose radii, the density radius and a shared density remain
unverified; EMPIAR-10076 is a heterogeneous assembly mixture. Failure to exclude
zero or to improve interval width is retained. Do not claim this additional
cohort resolves the first review's R1 or R2 by itself.
