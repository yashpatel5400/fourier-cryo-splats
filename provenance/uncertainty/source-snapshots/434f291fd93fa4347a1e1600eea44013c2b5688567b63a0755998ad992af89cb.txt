# Discrete-view feasibility results

30 September 2026. The repaired, unchanged design completes in 55.24 seconds.
All three stacks, sixteen repeats, four candidate maps and both noise settings
are retained. Of 384 mixture-envelope optimizations, 382 meet the 1e-4
log-likelihood-gap target; the largest remaining gap is 0.0005934. Every ratio
uses the upper bound, including nonconverged cases. The repaired first repeat's
means, observations and latent indices exactly reproduce the failed attempt.

Counts below are detections of a false candidate among 16 simulated repeats,
at e-value threshold 20. They are not precise power estimates or experimental
coverage measurements. The numerator is the oracle true discrete mixture.

| Stack | 50% local attenuation: known law / unknown shared law | Full local attenuation: known law / unknown shared law | Full attenuation: independent pose profiling |
|---|---:|---:|---:|
| 10028 | 15 / 0 | 16 / 14 | 0 |
| 10049 | 15 / 0 | 16 / 9 | 0 |
| 10076 | 5 / 0 | 16 / 0 | 0 |

These rows use the supplied noise variance. Full-region attenuation changes
the unit-norm densities by L2 norms 0.0863, 0.2438 and 0.1062, respectively;
mean squared whitened projection changes are 0.8363, 0.7848 and 0.1976.
Thus equal attenuation percentages are not equal statistical effects.
The zero-signal candidate is detected in every repeat by all comparators and
both noise treatments; that easy control alone is not evidence of utility.

Profiling a separate unrestricted positive variance inside every image/cell
envelope yields **zero local-attenuation detections by the unknown-shared-law
test on all three stacks**. At the true map, its median log-e-value penalties
are about 98–103, versus 34–40 with known variance. The fixed-known-law
comparator also loses most local detections. Its profiled-noise version is
an envelope comparator, not a normalized likelihood for a known variance.

The results separate two costs: assigning each image its own best pose is
especially permissive, while profiling one shared viewing distribution can
retain some structural discrimination. Relaxing variance independently at
every image/cell loses considerably more information. The matched true-map
ratios are at most one by construction with an oracle numerator; their zero
rejection counts are an algebra control, not calibration evidence.

**Decision.** There is limited evidence to investigate a common unknown noise
scale while retaining the same discrete catalog. This cannot fix the third
stack's known-noise deficiency by itself. No continuous-pose implementation
or experimental claim is justified yet. The new shared-scale calculation will
be post-outcome development on the existing simulated observations, not a
fresh confirmation or a selected replacement for these negative outcomes.
