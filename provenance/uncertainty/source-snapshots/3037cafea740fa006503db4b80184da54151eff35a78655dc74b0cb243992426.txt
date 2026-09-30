# Experimental application of all completed cubic estimators

30 September 2026. Declared after registered-target-sensitivity-v1, before
computing the three experimental intervals below. The registered generator
raised the original cubic estimator's minimum simulated sign power to 0.95725.
This motivates a development application, not a fresh confirmatory study.

Use EMPIAR-10049, the previously selected pilot_region_1 with Gaussian target
SD 20 A, and all three completed cubic-weight-probe, cubic-coordinate-probe,
and cubic-subspace-probe estimators. Do not refit weights, select a winner,
change target, or change their one-degree/0.5 A joint pose class, density
radius two, pilot norm one, enclosing domain, or audited bias bound. The
subspace estimator is numerically identical to the original; its tiny bias
difference is an audit difference, not a new fit improvement.

Use the existing inference_half0 geometry (128 particles, radius 12, seed
609315), original independent pilot amplitude, and all 128 previously
reserved noise-calibration-v1 exposure representatives. These pixels have
already been analyzed. Verify their published download hashes and disjoint
exposure identifiers before use. The spectral error allowance is 1e-6/12
per audit, the noise-calibration failure allowance is .005/12, and the
Gaussian tail allowance is (.045-1e-6)/12. Thus each conditional procedure
has nominal error allowance .05/12; union bounding these three fixed
procedures gives at most .0125 under all required independence/model
conditions. This is not a simultaneous guarantee with all historical
adaptively explored alternatives. Data reuse and unverified physical
conditions preclude claiming experimental coverage.

First replay the original fixed-pose estimator's raw observed center and its
published fresh-calibration SD and half-width to relative/absolute tolerance
1e-9. Then uncenter each cubic weight vector, apply the existing group-noise
variance upper bound, and form the bias-aware interval with its saved final
joint bias bound. Use the unchanged no-data fallback rule. Retain calibrated
SD, center, interval, zero exclusion, fallback, and both original-frame and
registered-reference inclusion. Both references use their published pilot-only
amplitude scaling. Reference inclusion is agreement, not coverage validation.
Also retain the published supplied-noise SD to expose the scale discrepancy.

Archive all inputs, source snapshots, calibration arrays, observations, and
the three raw weight arrays. Do not overwrite existing records. A 15-minute
wall budget applies to the complete application. Preserve failures. No new
particle download, optimizer run, or inference-based registration is allowed.
