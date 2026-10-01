# Covariance-aware scores and equal-budget calibration — development checkpoint

The 33-page ICML-format paper adds a fresh classical covariance-aware score comparison
and a diagnostic of noise-replica allocation on all three fitted Gaussian maps.
It retains failures and every correct-null control. All three full Claude Fable
5.1 reviews remain rejections; there is no new full acceptance assessment.

Fisher scoring improves the largest 10028 combined-moment case, with a
conditional projected-power difference interval excluding zero, but smaller
changes and 10076 remain difficult. The result does not establish a general
power advantage or experimentally calibrated regional-density uncertainty.

The allocation comparison holds conditional noise-draw count fixed while
changing views/replicas, and evaluates both allocations on fresh common test
images. Estimated envelopes decrease, but increased confidence slack often
prevents a power gain. It is not an optimized allocation, an equal-wall-time comparison or
a repeated-calibration error experiment. No unadjusted minimum across methods,
scores or allocations is used as a valid combined test.

All arrays are supplied in four numerical bundles, split by stack to keep each
asset under 2GB. These include the full noisy training inputs, not only summary
statistics. Twenty-five targeted tests pass. Separate implementations verify
covariances, constrained score directions, probability bounds, critical values,
all projections/group vectors and 9,450 conservative difference intervals.

See `research/uncertainty/REPRODUCE-REVISION9.md` for frozen commits, commands,
checks and dependencies. Existing releases and unaltered reviews remain public.
