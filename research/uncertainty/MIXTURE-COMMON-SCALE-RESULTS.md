# Shared-scale discrete likelihood: complete development follow-up

30 September 2026. Protocol and implementation were published in `fdd4cd7`
before this run. This reuses every simulated image from the preceding screen;
it is post-outcome method development, not independent confirmation. The model
still has only 64 allowed views and an oracle normalized numerator. It does
not validate continuous orientations or an experimental noise model.

All 192 cases completed in 72.826 seconds under concurrent Mac workloads.
Every outer likelihood bracket met the declared one-log-unit tolerance; the
largest gap was 0.936776, with at most three scale-interval nodes per case.
This is ordinary floating-point evidence, not certified interval arithmetic.
The exact endpoint mixture fits, anchors, weights and stopping statuses are
retained in the per-stack JSON. Outer convergence does not assert that every
inner mixture fit met its tighter iteration tolerance.

The test rejects at a log ratio of at least log(20). Counts below are out of
16 paired simulations per candidate and stack, not 16 independent biological
datasets. All methods use the same unknown shared viewing-law family.

| Stack | Full local removal, supplied scale | Full removal, shared unknown scale | Half local removal, shared scale | Zero signal, shared scale |
|---|---:|---:|---:|---:|
| EMPIAR-10028 | 14 | 14 | 0 | 16 |
| EMPIAR-10049 | 9 | 8 | 0 | 16 |
| EMPIAR-10076 | 0 | 0 | 0 | 16 |

The median log ratios for full local removal with the shared scale are
20.8912, 3.30983 and -28.1554. The corresponding true-map medians are
-40.2265, -35.3235 and -34.7782. True-map nonrejection here is a pointwise
algebraic consequence of the matched oracle numerator and enlarged null;
it is not an empirical calibration result. Independent scale profiling in
each image/cell had removed all local detections in the earlier relaxation.
Retaining a shared scale therefore matters in this restricted calculation.

Relative to the supplied-scale calculation, median reductions in log ratio
across nonzero candidate maps are 0.172994, 0.599882 and 0.513020; maxima are
3.88118, 4.83541 and 1.19821. These include numerical bracket slack. All cases
passed the check that the shared-scale upper bound exceeds a feasible
supplied-scale likelihood, with the recorded numerical gaps accounted for.

The failed initial loose-bound numerical test remains archived separately.
No observation, candidate, repeat or stack was removed. The source arrays and
their hashes are those of `mixture-validation-preflight-v2`; no new simulation
draws were generated. Results are in `mixture-common-scale-v1`.

**Next scientific gate.** A practical numerator learned independently of the
validation sample and an upper bound covering continuous orientations are
both missing. More oracle finite-catalog repetitions alone would not resolve
either requirement. Shared scale also does not cover colored noise, incorrect
CTFs, image dependence or molecular heterogeneity. The pooled experimental
EMPIAR-10076 stack cannot acquire a homogeneous interpretation from this
single-reference simulation.
