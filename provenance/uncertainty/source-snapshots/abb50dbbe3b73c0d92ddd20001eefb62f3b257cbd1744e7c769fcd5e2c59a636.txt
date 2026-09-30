# Registered-reference target sensitivity with unchanged estimators

30 September 2026. Declared after the ten reference-registration FSC outcomes.
This is post-outcome sensitivity to the reference frame, not new confirmation,
calibration or an estimator improvement. Use all three published transforms
from reference-registration-v1 without any further alignment or parameter fit.
Preserve every original reference-generator and experimental comparison.

For all twelve locked 20 A features, retain the published fixed-pose weights,
the fixed/shift-only/one-degree/two-degree enclosing-domain half-widths,
fallback rules, supplied simulation noise and family budgets. Normalize the
registered 64-cube reference to unit L2 constant-cell norm, exactly as for the
old simulation generator. Recompute its target and expected affine center at
the nominal, coherent-x and seeded random-boundary configurations used by the
old pose study (seed 610281+dataset; the fixed-pose case needs only nominal).
Retain every coverage/sign-power value as an analytic model calculation.
Before the changed-frame calculation, replay each original nominal reference
center against its saved value; reject a mismatch rather than attributing it
to registration. Reuse geometry count 128, radius 12 and seed 609315. No weights,
targets, uncertainty bounds, noise calibration, or observed intervals change.

On 10049's already selected region 1, repeat the same three configurations for
each completed optimized cubic weight vector: original, coordinate and subspace.
Use each saved half-width/no-data fallback and noise SD. These three alternative
estimators are reported separately; do not select the most favorable one or
claim simultaneous validity across that post-outcome collection.

For the existing 48 frozen fresh-noise experimental intervals, compare the
registered reference after fitting only its amplitude to the same old 256
pilot particles at radius 5 and seed 609681 used by the original amplitude
calculation. The pilot-density amplitude and interval centers/widths stay fixed.
Replay the original pilot-density amplitude as a consistency check. Report
the newly fitted reference amplitude, its sign, the reference-to-pilot L2
distance and whether it lies in the supplied radius-two class; do not force a
negative amplitude or out-of-class reference into the model. All 48 original
unregistered target values/inclusion flags remain alongside the new values.
Pilot amplitude fitting is a descriptive comparison convention, not an estimate
of true experimental density or a calibration guarantee.

Run every declared feature/pose/estimator case, retaining failures and no-data
fallbacks. Save signal vectors, unit-norm generators, geometry, indices, targets,
all scalar outcomes and exact source/input hashes. There is no particle refit,
new noise sample, threshold selection or new target selection. Budget 1,800
seconds per geometry, checked between signal evaluations; keep any timed-out
geometry and continue the remaining ones. Use single-thread FINUFFT on this Mac.
