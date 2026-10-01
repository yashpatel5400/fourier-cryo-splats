# Experimental application of enriched and joint cubic estimators

30 September 2026. Declared after completion of the small adaptive enrichment
and during the first joint-design fit, before either new estimator is applied
to the experimental pixels. Earlier fixed/three-cubic applications and their
raw/centered/projected noise outcomes have been seen. This is development on
reused data, not a new frozen confirmation experiment.

Apply both completed weight vectors, with no refitting or selection by image
outcome, to the unchanged 10049 region-1 target: 20 A SD, 128 radius-12 inference
particles, one-degree/0.5 A joint pose class, B=2 and pilot norm one. Require
completed successful fit records and verify the exact arrays, indices, scales,
noise and pilot/target metadata before use. A missing/failed fit is not a
successful estimator and cannot be silently replaced by another checkpoint.

Use the same 128 exposure-distinct calibration representatives and original
128 inference particles (26 exposures, maximum group size twelve). Replay
the frozen original fixed-estimator center/SD/half-width and all three older
cubic raw/centered SD bounds before any new outcome. Use the unchanged grouped
noncentral Gaussian trace bound and its fixed 127-row Helmert version as two
separately reported procedures. The inference variance bound permits arbitrary
within-exposure Gaussian dependence by inflating each weight row by sqrt(group
size). Common covariance, noise/metadata independence, pose radii, support and
density radius remain unverified experimental assumptions.

For each of the four new estimator/calibration combinations, keep the same
alpha_noise=(.045-1e-6)/12, beta=.005/12, delta=1e-6/12. Recompute the folded-
normal interval with its unchanged continuous bias. If it exceeds no data,
recenter at the pilot and use the no-data radius. Retain both original and
pilot-registered approximate-reference comparisons; neither is a coverage
label. Do not take an unadjusted minimum across these procedures, select a
favorable calibration, or infer simultaneous experimental validity from a
hypothetical fixed-design union bound after this adaptive development.

Report all four centers, SD bounds, bias terms, widths, fallbacks, zero
exclusions and reference disagreements. Preserve saved transformed calibration
and inference observations and raw-frame weights for replay; disclose those
observations in any artifact release. Use a 900-second application boundary
budget, no model fitting and no additional data collection. The only outcome
selection is the previously fixed estimator's own completed design/audit rule.
