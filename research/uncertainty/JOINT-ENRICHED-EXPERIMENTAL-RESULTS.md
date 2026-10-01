# Completed enriched and joint experimental application

30 September 2026. Both new estimators are applied under the committed
four-interval protocol to the same previously inspected 10049 images. The
runner first replays the frozen fixed-estimator center, SD upper and width,
and the raw/centered SDs for all three earlier cubic designs. It does not fit
new weights to the calibration or inference pixels.

| Design | Calibration | Center | Half-width | SD upper |
|---|---|---:|---:|---:|
| Enriched | Raw | 2.463500 | 4.282612 | 1.009318 |
| Enriched | Centered | 2.463500 | 3.060135 | .552112 |
| Joint | Raw | 2.468957 | 4.219517 | .997224 |
| Joint | Centered | 2.468957 | 3.010881 | .545193 |

All four intervals contain zero, the registered approximate reference 2.488459
and the original-frame reference .710119. None uses the no-data fallback.
Thus the modest simulation-width improvements do not resolve the experimental
sign-detection limitation. These SD upper bounds include signal energy and
conservatism; they are not measurements of pure noise.

There are 128 calibration exposure representatives; the 128 inference particles
occupy 26 exposures with maximum group size twelve. Each procedure retains
noise, calibration and spectral allocations (.045 - 1e-6)/12, .005/12 and
1e-6/12. Raw and centered alternatives are reported separately, without taking
an unadjusted minimum. The data have been inspected previously; this is not
fresh experimental validation or a simultaneous-coverage claim across the
development history. Common Gaussian covariance, independent exposures and
metadata, and the density/pose radii remain unverified assumptions.

The computation takes .383 seconds after fitting and reuses the recorded arrays.
`joint-enriched-experimental-v1/10049.json` retains every interval and replay;
its integrity record verifies all 23 inputs, 70 source snapshots and the output
arrays. Figure `cubic-design-comparison` displays all five cubic designs, both
simulation frames, and both experimental calibration procedures.

